use std::fs::File;
use std::io::{BufReader, Read, BufRead};

use serde::{Serialize, Deserialize};

use rand::Rng;
use rand::prelude::*;

use crate::io;

// Not as well implemented, but based off of ELECTRA model
// https://github.com/google-research/electra
pub struct DiscriminatorMasked {
    pub kmers: Vec<Vec<u8>>,
    pub id: String,
    pub taxons: Vec<usize>,
    pub taxon: usize,
    pub truth: Vec<u8>,
}

pub struct DiscriminatorMaskedGenerator {
    kmer_window_generator: KmerWindowGenerator,
    replacement_pct: f32,
    rng: rand::prelude::ThreadRng,
    k: usize,
    // TODO: Put alphabet here so we can train proteins as well...
}

pub struct KmerWindow {
    pub kmers: Vec<Vec<u8>>,
    pub id: String,
    pub taxons: Vec<usize>,
    pub taxon: usize,
}

pub struct KmerWindowGenerator {
    sequences: io::Sequences,
    window_size: usize,
    k: usize,
    kmer_generator: Kmers,
    needed_sequence: usize,
    curseq: io::Sequence,
    offset: usize,
}

impl DiscriminatorMaskedGenerator {
    pub fn new(
        replacement_pct: f32, 
        k: usize,
        generator: KmerWindowGenerator) -> DiscriminatorMaskedGenerator 
    {
        let mut rng = rand::thread_rng();

        DiscriminatorMaskedGenerator {
            kmer_window_generator: generator,
            replacement_pct,
            rng,
            k,
        }
    }
}

impl Iterator for DiscriminatorMaskedGenerator {
    type Item = DiscriminatorMasked;

    fn next(&mut self) -> Option<DiscriminatorMasked> {
        let mut next_item = match self.kmer_window_generator.next() {
            Some(x) => x,
            None    => return None
        };

        let KmerWindow { mut kmers, id, taxons, taxon } = next_item;

        // TODO: Make switchable, so we can train protein sequences
        // ~2% chance of an N
        let choices = [(b'A', 48), (b'C', 48), (b'T', 48), (b'G', 48), (b'N', 4)]; 

        let mut truth: Vec<u8> = Vec::with_capacity(kmers.len());

        for kmer in kmers.iter_mut() {
            if self.rng.gen::<f32>() < self.replacement_pct {
                let mut new_kmer: Vec<u8> = Vec::with_capacity(self.k);
                for _i in 0..self.k {
                    new_kmer.push(choices.choose_weighted(&mut self.rng, |item| item.1).unwrap().0);
                }
                *kmer = new_kmer;
                truth.push(1);
            } else {
                truth.push(0);
            }
        }

        return Some(DiscriminatorMasked { kmers, id, taxons, taxon, truth })
    }
}

impl KmerWindowGenerator {
    // TODO: If RC, then do RC here.... or maybe elsewhere?
    // TODO: Make offset automatic? Or a diff fn to handle it internally?

    pub fn new(filename: String,
                k: usize,
                window_size: usize,
                offset: usize,
            ) -> KmerWindowGenerator {
        let mut sequences = io::Sequences::new(filename);
        let curseq = match sequences.next() {
            Some(x) => x,
            None    => panic!("File is empty or invalid format!")
        };
        
        let kmer_generator = Kmers::new(curseq.seq.clone(), k, offset);
        let needed_sequence = k * window_size;

        KmerWindowGenerator {
            sequences,
            window_size,
            k,
            kmer_generator,
            needed_sequence, 
            curseq,
            offset,
        }
    }
}

// TODO: Set up a function to apply to the output of KmerWindowGenerator 
// that will randomly mask, or randomly replace some percent

impl Iterator for KmerWindowGenerator {
    type Item = KmerWindow;

    fn next(&mut self) -> Option<KmerWindow> {
        // While instead of if, because if we get a too short sequence we should skip it...
        while (self.kmer_generator.len - self.kmer_generator.curpos) < self.needed_sequence {
            // Don't have enough remaining sequence to fill out this.
            // Get the next sequence...            

            let curseq: io::Sequence = match self.sequences.next() {
                Some(x) => x,
                None    => return None // That's it... no more!
            };
            self.curseq = curseq;
            self.kmer_generator.seq = self.curseq.seq.clone();
            self.kmer_generator.curpos = 0;
            self.kmer_generator.len = self.kmer_generator.seq.len();
        }

        let mut kmers: Vec<Vec<u8>> = Vec::with_capacity(self.window_size);
        for _i in 0..self.window_size {
            kmers.push(self.kmer_generator.next().unwrap());
        }

        Some(KmerWindow { 
            kmers, 
            id: self.curseq.id.clone(), 
            taxons: self.curseq.taxons.clone(),
            taxon: self.curseq.taxon.clone(),
        })

    }
}

pub struct Kmers {
    seq: Vec<u8>,
    k: usize,
    offset: usize, 
    curpos: usize,
    len: usize,
}

impl Kmers {
    pub fn new(
        seq: Vec<u8>,
        k: usize, 
        offset: usize,
    ) -> Kmers {
        let len = seq.len();
       
        Kmers { 
            seq, k, offset, len, curpos: 0
        }
    }
}

impl Iterator for Kmers {
    type Item = Vec<u8>;

    #[inline(always)]
    fn next(&mut self) -> Option<Vec<u8>> {
        if (self.offset + self.curpos + self.k) > self.len {
            None
        } else {
            let start = self.offset + self.curpos;
            let end = self.offset + self.curpos + self.k;
            self.curpos = self.curpos + self.k;
            Some(self.seq[start..end].to_vec())
        }
    }
}