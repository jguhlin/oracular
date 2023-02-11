// Std crates
use std::collections::HashMap;
use std::f32::consts::E;
use std::sync::atomic::AtomicBool;
use std::sync::atomic::Ordering;
use std::sync::Arc;
use std::thread;
use std::thread::{park, JoinHandle};
use std::io::Write;

// External crates
use crossbeam::queue::ArrayQueue;
use crossbeam::utils::Backoff;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rand::prelude::*;
use rand_xoshiro::rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;

// So this can be used from python without running the database twice
pub use acc2tax::*;

// Internal Crates
use super::*;
use liboracular::io;
use liboracular::kmers::KmerWindow;
use liboracular::kmers::{
    rc_kmerwindow, replace_blank, replace_random, DiscriminatorMasked,
    DiscriminatorMaskedGenerator, Gff3Kmers, Gff3KmersIter, KmerCoordsWindow, KmerCoordsWindowIter,
    KmerWindowGenerator,
};
use libsfasta::prelude::*;

#[pyclass(frozen)]
pub struct TaxonomyReturn {
    #[pyo3(get)]
    pub sequence: Vec<Vec<bool>>,
    #[pyo3(get)]
    pub taxonomy: Vec<Vec<bool>>,
}

#[derive(Clone)]
#[pyclass(frozen)]
pub struct Taxon {
    #[pyo3(get)]
    pub id: usize,
    #[pyo3(get)]
    pub name: String,
    #[pyo3(get)]
    pub seq_count: usize,
    #[pyo3(get)]
    pub seq_length: usize,
}

// ** Taxonomic Classification **
#[pyclass]
pub struct TaxonomicSequenceGenerator {
    // queueimpl: QueueImpl<TaxonomyReturn>,
    #[pyo3(get)]
    pub child_taxons: Vec<Taxon>,
}

impl Drop for TaxonomicSequenceGenerator {
    fn drop(&mut self) {
        // self.queueimpl.shutdown();
    }
}

#[pymethods]
impl TaxonomicSequenceGenerator {
    #[new]
    fn new(
        k: usize,
        filename: String,
        window_size_min: usize,
        window_size_max: usize,
        taxon_root: usize,
        queue_size: usize,
        threads: usize,
        seed: usize,
        acc2tax_threads: usize,
        acc2tax_filename: String,
        acc2tax_nodes_filename: String,
        acc2tax_names_filename: String,
    ) -> PyResult<Self> {

        assert!(window_size_min <= window_size_max);
        assert!(window_size_min > 0);
        assert!(window_size_max > 0);

        acc2tax::init(
            acc2tax_threads,
            acc2tax_filename,
            acc2tax_nodes_filename,
            acc2tax_names_filename,
        );

        let ranks = acc2tax::get_taxon_ranks();
        // if ranks.contains(&taxon_rank) {
        // return Err(PyValueError::new_err(format!("Taxon rank {} not found",
        // taxon_rank))) }

        let child_taxons = acc2tax::get_child_taxons_names(taxon_root);
        if child_taxons.is_empty() {
            return Err(PyValueError::new_err(format!(
                "Taxon root {} has {} children",
                taxon_root,
                child_taxons.len()
            )));
        }

        let child_taxon_rank = acc2tax::get_taxon_rank(child_taxons[0].0);

        let mut balancer: HashMap<usize, usize> =
            child_taxons.iter().map(|(taxon, _)| (*taxon, 0)).collect();
        let mut totals: HashMap<usize, usize> =
            child_taxons.iter().map(|(taxon, _)| (*taxon, 0)).collect();
        let mut total_seq_len: HashMap<usize, usize> =
            child_taxons.iter().map(|(taxon, _)| (*taxon, 0)).collect();

        assert!(!child_taxons.is_empty());
        println!("Identified {} child taxons", child_taxons.len());

        let mut buf = std::io::BufReader::new(std::fs::File::open(filename.clone()).unwrap());
        let mut sfasta = SfastaParser::open_from_buffer(&mut buf, true).unwrap();

        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed as u64);

        // This preloads the seqlocs if not already loaded
        // sfasta.get_seqlocs().expect("Unable to get seqlocs");

        println!("Identifying suitable sequences...");
        let min_length = window_size_min * k;

        let idx = sfasta.seqlocs.as_ref().unwrap().index.as_ref().unwrap();
        let valid_indices = idx.iter().enumerate().filter_map(|(i, seqloc)| 
            { 
                let seqlen = sfasta.seqloc_len_loaded(seqloc);
                if seqlen >= min_length {
                    Some((i, seqloc, seqlen))
                } else {
                    None
                }
            }
        )
        .map(|(i, seqloc, seqlen)| {
            let id = sfasta.get_id_loaded(&seqloc.ids.as_ref().unwrap()).expect("Unable to get id Locs");
            if Python::with_gil(|py| py.check_signals()).is_err() {
                panic!("Interrupted");
            };
            (i, seqloc, seqlen, id)
        })  
        .filter_map(|(i, _seqloc, seqlen, id)| 
            {
            if Python::with_gil(|py| py.check_signals()).is_err() {
                panic!("Interrupted");
            };
            // Get taxons and check if any match the children
            let taxon = acc2tax::get_taxon(id.clone());
            let taxons = acc2tax::get_complete_taxonomy(taxon as usize);
            if let Some(x) = taxons.iter().filter(|x| balancer.contains_key(x)).next() {
                return Some((*x, i, seqlen, id))
            }
            None
            }
        ).collect::<Vec<(usize, usize, usize,String)>>();

        println!("Identified {} valid sequences", valid_indices.len());
        valid_indices.iter().for_each(|(taxon, _i, seqlen, _id)| {
            totals.entry(*taxon).and_modify(|x| *x += 1);
            total_seq_len.entry(*taxon).and_modify(|x| *x += seqlen);
        });

        // let taxon_indices: HashMap<usize, Vec<usize>> = HashMap::new();
        let taxon_indices = valid_indices
            .iter()
            .fold(HashMap::new(), |mut acc, (taxon, i, _seqlen, _id)| {
                acc.entry(*taxon).or_insert_with(Vec::new).push(*i);
                acc
            });

        println!("{:#?}", totals);

        let child_taxons: Vec<Taxon> = child_taxons
            .into_iter()
            .map(|(taxon, name)| {
                let seq_count = totals.get(&taxon).unwrap();
                let seqlen = total_seq_len.get(&taxon).unwrap();
                Taxon {
                    id: taxon,
                    name,
                    seq_count: *seq_count,
                    seq_length: *seqlen,
                }
            })
            .collect();

        // Write the taxon counts to a file
        let mut f = std::fs::File::create("taxon_counts.csv").unwrap();
        for taxon in child_taxons.iter() {
            writeln!(f, "{}\t{}\t{}\t{}", taxon.id, taxon.name, taxon.seq_count, taxon.seq_length).unwrap();
        }

        Ok(Self { child_taxons })
    }
}
