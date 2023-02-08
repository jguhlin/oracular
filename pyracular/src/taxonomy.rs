// Std crates
use std::sync::atomic::Ordering;
use std::thread;
use std::thread::{park, JoinHandle};
use std::sync::atomic::AtomicBool;
use std::sync::Arc;
use std::collections::HashMap;

// External crates
use crossbeam::queue::ArrayQueue;
use crossbeam::utils::Backoff;
use rand::prelude::*;
use rand_xoshiro::rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;

// So this can be used from python without running the database twice
pub use acc2tax::*;

// Internal Crates
use liboracular::io;
use liboracular::kmers::KmerWindow;
use liboracular::kmers::{
    rc_kmerwindow, replace_blank, replace_random, DiscriminatorMasked,
    DiscriminatorMaskedGenerator, Gff3Kmers, Gff3KmersIter, KmerCoordsWindow, KmerCoordsWindowIter,
    KmerWindowGenerator,
};
use libsfasta::prelude::*;
use super::*;

#[pyclass(frozen)]
pub struct TaxonomyReturn {
    #[pyo3(get)]
    pub sequence: Vec<Vec<bool>>,
    #[pyo3(get)]
    pub taxonomy: Vec<Vec<bool>>,
}

// ** Taxonomic Classification **
#[pyclass]
struct TaxonomicSequenceGenerator {
    // queueimpl: QueueImpl<TaxonomyReturn>,
    child_taxons: Vec<(usize, String)>,
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
        acc2tax::init(acc2tax_threads, acc2tax_filename, acc2tax_nodes_filename, acc2tax_names_filename);

        let ranks = acc2tax::get_taxon_ranks();
        // if ranks.contains(&taxon_rank) {
            // return Err(PyValueError::new_err(format!("Taxon rank {} not found", taxon_rank)))
        // }

        let child_taxons = acc2tax::get_child_taxons_names(taxon_root);
        if child_taxons.is_empty() {
            return Err(PyValueError::new_err(format!("Taxon root {} has {} children", taxon_root, child_taxons.len())))
        }

        let child_taxon_rank = acc2tax::get_taxon_rank(child_taxons[0].0);

        let mut balancer: HashMap<usize, usize> = child_taxons.iter().map(|(taxon, _)| (*taxon, 0)).collect();

        let mut buf = std::io::BufReader::new(std::fs::File::open(filename.clone()).unwrap());
        let mut sfasta = SfastaParser::open_from_buffer(&mut buf, false).unwrap();

        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed as u64);
        
        // This preloads the seqlocs if not already loaded
        sfasta.get_seqlocs().expect("Unable to get seqlocs");

        let mut valid_indices = Vec::new();

        let min_length = window_size_min * k;

        for i in 0..sfasta.seqlocs.as_ref().unwrap().total_seqlocs {
            let seqloc = sfasta.get_seqloc(i).expect("Unable to get seqloc");
            if let Some(seqloc) = seqloc {
                // If the sequence is too small we skip it
                if sfasta.seqloc_len(&seqloc) < min_length {
                    continue;
                }

                // Sequence must also be part of a child taxon
                if let Some(ids_loc) = &seqloc.ids {
                    let id = sfasta.get_id(ids_loc).expect("Unable to get id");
                    // Get taxons and check if any match the children
                    let taxon = acc2tax::get_taxon(id);
                    let taxons = acc2tax::get_complete_taxonomy(taxon as usize);
                    if taxons.iter().any(|x| balancer.contains_key(x)) {
                        valid_indices.push(i);
                    }
                }
            }
        }

        Ok(Self {
            child_taxons,
        })

    }
}

