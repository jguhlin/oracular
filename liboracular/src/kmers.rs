use std::fs::File;
use std::io::{BufReader, Read, BufRead};

use serde::{Serialize, Deserialize};

use rand::Rng;
use rand::prelude::*;

use crate::io;
use crate::sfasta;
use crate::utils;

// Not as well implemented, but based off of ELECTRA model
// https://github.com/google-research/electra
#[derive(Debug, PartialEq)]
pub struct DiscriminatorMasked {
    pub kmers: Vec<Vec<u8>>,
    pub id: String,
    pub truth: Vec<u8>,
}

/*pub struct DiscriminatorMaskedA2T {
    pub kmers: Vec<Vec<u8>>,
    pub id: String,
    pub taxons: Vec<usize>,
    pub taxon: usize,
    pub truth: Vec<u8>,
} */

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
//    pub taxons: Vec<usize>,
//    pub taxon: usize,
}

pub struct KmerWindowGenerator {
    sequences: Box<dyn Iterator<Item = io::Sequence>>,
    window_size: usize,
    k: usize,
    kmer_generator: Kmers,
    needed_sequence: usize,
    curseq: io::Sequence,
    offset: usize,
    rc: bool,
}

impl DiscriminatorMaskedGenerator {
    pub fn new(
        replacement_pct: f32, 
        k: usize,
        generator: KmerWindowGenerator) -> DiscriminatorMaskedGenerator 
    {
        let rng = rand::thread_rng();

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
        let next_item = match self.kmer_window_generator.next() {
            Some(x) => x,
            None    => return None
        };

        let KmerWindow { mut kmers, id } = next_item;

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

        // return Some(DiscriminatorMasked { kmers, id, taxons, taxon, truth })
        return Some(DiscriminatorMasked { kmers, id, truth })
    }
}

impl KmerWindowGenerator {

    pub fn new(filename: String,
                k: usize,
                window_size: usize,
                offset: usize,
                rc: bool,
            ) -> KmerWindowGenerator {
        let mut sequences = Box::new(sfasta::Sequences::new(filename));
        let mut curseq = match sequences.next() {
            Some(x) => x,
            None    => panic!("File is empty or invalid format!")
        };

        if rc {
            let io::Sequence { id, mut seq } = curseq;

            utils::complement_nucleotides(&mut seq);
            seq.reverse();

            curseq = io::Sequence { 
                id, 
                seq,
            };
        }
        
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
            rc,
        }
    }
}

// TODO: Set up a function to apply to the output of KmerWindowGenerator 
// that will randomly mask, or randomly replace some percent

impl Iterator for KmerWindowGenerator {
    type Item = KmerWindow;

    fn next(&mut self) -> Option<KmerWindow> {
        // While instead of if, because if we get a too short sequence we should skip it...
        while (self.kmer_generator.len - self.kmer_generator.curpos) <= self.needed_sequence {
            // Don't have enough remaining sequence to fill out this.
            // Get the next sequence...            

//            println!("Getting new sequence...");
            let mut curseq: io::Sequence = match self.sequences.next() {
                Some(x) => x,
                None    => { println!("Reporting we are out of sequence!"); return None }  // That's it... no more!
            };
            
//            println!("Getting new sequence... {}", &curseq.id);

            if self.rc {
                let io::Sequence { id, mut seq } = curseq;
    
                utils::complement_nucleotides(&mut seq);
                seq.reverse();
    
                curseq = io::Sequence { 
                    id, 
                    seq,
                };
            }

            self.curseq = curseq.clone();
            let kmer_generator = Kmers::new(curseq.seq.clone(), self.k, self.offset);
	    self.kmer_generator = kmer_generator;

//            self.kmer_generator.seq = self.curseq.seq.clone();
//            self.kmer_generator.curpos = 0;
//            self.kmer_generator.len = self.kmer_generator.seq.len();
        }

        let mut kmers: Vec<Vec<u8>> = Vec::with_capacity(self.window_size);
        for _ in 0..self.window_size {
            match self.kmer_generator.next() {
                Some(x) => {
                    if x.len() == self.k {
                        kmers.push(x)
                    } else {
                        panic!("Invalid KMER");

                    }
                },
                None    => return None
            };
        }

        Some(KmerWindow { 
            kmers, 
            id: self.curseq.id.clone(), 
//            taxons: self.curseq.taxons.clone(),
//            taxon: self.curseq.taxon.clone(),
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
        if (self.offset + self.curpos + self.k) >= self.len {
            None
        } else {
            let start = self.offset + self.curpos;
            let end = self.offset + self.curpos + self.k;
            self.curpos = self.curpos + self.k;
            Some(self.seq[start..end].to_vec())
        }
    }
}
