use rand::Rng;
use rand::prelude::*;
use std::fmt;

use crate::io;
use crate::sfasta;
use crate::utils;
use crate::gff3;
use crate::intervals;

type Kmer = Vec<u8>;
type Coords = (usize, usize);

pub struct KmerCoords {
    pub kmer: Vec<u8>,
    pub coords: Coords,
}

impl fmt::Display for KmerCoords {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {}-{}", 
            std::str::from_utf8(&self.kmer).expect("Unable to convert Kmer to Str"),
            self.coords.0,
            self.coords.1,
        )
    }
}

impl fmt::Debug for KmerCoords {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, 
            "{}: {}-{}", 
            std::str::from_utf8(&self.kmer).expect("Unable to convert Kmer to Str"),
            self.coords.0,
            self.coords.1,
        )
    }
} 

impl std::cmp::PartialEq<&[u8]> for KmerCoords {
    fn eq (&self, other: &&[u8]) -> bool {
        self.kmer == *other
    }
}

#[derive(Debug)]
pub struct KmerWindow {
    pub kmers: Vec<Kmer>,
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

// Loosely based off of ELECTRA model
// https://github.com/google-research/electra
#[derive(Debug, PartialEq)]
pub struct DiscriminatorMasked {
    pub kmers: Vec<Vec<u8>>,
    pub id: String,
    pub truth: Vec<u8>, // TODO: Should this be bool? 
                        // Maybe it's not because doesn't work with python? 
                        // Or because tf works best with 0 and 1?
}

pub struct DiscriminatorMaskedGenerator {
    kmer_window_generator: KmerWindowGenerator,
    replacement_pct: f32,
    rng: rand::prelude::ThreadRng,
    k: usize,
    // TODO: Put alphabet here so we can train proteins as well...
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
        Some(DiscriminatorMasked { kmers, id, truth })
    }
}

impl KmerWindowGenerator {
    pub fn new(filename: String,
                k: usize,
                window_size: usize,
                offset: usize,
                rc: bool,
            ) -> KmerWindowGenerator {
        let sequences = Box::new(sfasta::Sequences::new(filename));
        let mut sequences = Box::new(io::SequenceSplitter3N::new(sequences));
        let mut curseq = match sequences.next() {
            Some(x) => x,
            None    => panic!("File is empty or invalid format!")
        };

        if rc {
            let io::Sequence { id, mut seq, location } = curseq;

            utils::complement_nucleotides(&mut seq);
            seq.reverse();

            // TODO: How to deal with rc?

            curseq = io::Sequence { 
                id, 
                seq,
                location,
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

impl Iterator for KmerWindowGenerator {
    type Item = KmerWindow;

    fn next(&mut self) -> Option<KmerWindow> {
        // While instead of if, because if we get a too short sequence we should skip it...
        while (self.kmer_generator.len - self.kmer_generator.curpos) <= self.needed_sequence {
            let mut curseq: io::Sequence = match self.sequences.next() {
                Some(x) => x,
                None    => return None  // That's it... no more!
            };
            
            if self.rc {
                let io::Sequence { id, mut seq, location } = curseq;
    
                utils::complement_nucleotides(&mut seq);
                seq.reverse();
    
                curseq = io::Sequence { 
                    id, 
                    seq,
                    location,
                };
            }

            self.curseq = curseq.clone();
            let kmer_generator = Kmers::new(curseq.seq.clone(), self.k, self.offset);
	        self.kmer_generator = kmer_generator;
        }

        let mut kmers: Vec<Vec<u8>> = Vec::with_capacity(self.window_size);
        let mut coords: Vec<Coords> = Vec::with_capacity(self.window_size);
        for _ in 0..self.window_size {
            match self.kmer_generator.next() {
                Some((x, c)) => {
                    if x.len() == self.k {
                        kmers.push(x);
                        coords.push(c);
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
        })
    }
}

#[derive(Debug)]
pub struct KmerCoordsWindow {
    pub kmers: Vec<Kmer>,
    pub coords: Vec<Coords>,
    pub id: String,
}

pub struct KmerCoordsWindowGenerator {
    sequences: Box<dyn Iterator<Item = io::Sequence>>,
    window_size: usize,
    k: usize,
    kmer_generator: Kmers,
    needed_sequence: usize,
    curseq: io::Sequence,
    offset: usize,
    rc: bool,
}

impl KmerCoordsWindowGenerator {
    pub fn new(filename: String,
                k: usize,
                window_size: usize,
                offset: usize,
                rc: bool,
            ) -> KmerCoordsWindowGenerator {
        let sequences = Box::new(sfasta::Sequences::new(filename));
        let mut sequences = Box::new(io::SequenceSplitter3N::new(sequences));
        let mut curseq = match sequences.next() {
            Some(x) => x,
            None    => panic!("File is empty or invalid format!")
        };

        if rc {
            let io::Sequence { id, mut seq, location } = curseq;

            utils::complement_nucleotides(&mut seq);
            seq.reverse();

            // TODO: How to deal with rc?

            curseq = io::Sequence { 
                id, 
                seq,
                location,
            };
        }
        
        let kmer_generator = Kmers::new(curseq.seq.clone(), k, offset);
        let needed_sequence = k * window_size;

        KmerCoordsWindowGenerator {
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

impl Iterator for KmerCoordsWindowGenerator {
    type Item = KmerCoordsWindow;

    fn next(&mut self) -> Option<KmerCoordsWindow> {
        // While instead of if, because if we get a too short sequence we should skip it...
        while (self.kmer_generator.len - self.kmer_generator.curpos) <= self.needed_sequence {
            let mut curseq: io::Sequence = match self.sequences.next() {
                Some(x) => x,
                None    => return None  // That's it... no more!
            };
            
            if self.rc {
                let io::Sequence { id, mut seq, location } = curseq;
    
                utils::complement_nucleotides(&mut seq);
                seq.reverse();
    
                curseq = io::Sequence { 
                    id, 
                    seq,
                    location,
                };
            }

            self.curseq = curseq.clone();
            let kmer_generator = Kmers::new(curseq.seq.clone(), self.k, self.offset);
	        self.kmer_generator = kmer_generator;
        }

        let mut kmers: Vec<Vec<u8>> = Vec::with_capacity(self.window_size);
        let mut coords: Vec<Coords> = Vec::with_capacity(self.window_size);
        for _ in 0..self.window_size {
            match self.kmer_generator.next() {
                Some((x, c)) => {
                    if x.len() == self.k {
                        kmers.push(x);
                        coords.push(c);
                    } else {
                        panic!("Invalid KMER");
                    }
                },
                None    => return None
            };
        }

        // TODO: Math is probably different when RC
        let coords = coords.iter().map(
                        |(x, y)| 
                        (x+self.curseq.location,
                         y+self.curseq.location,)).collect();

        Some(KmerCoordsWindow { 
            kmers, 
            id: self.curseq.id.clone(), 
            coords,
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
    type Item = (Vec<u8>, Coords);

    #[inline]
    fn next(&mut self) -> Option<(Vec<u8>, Coords)> {
        if (self.offset + self.curpos + self.k) >= self.len {
            None
        } else {
            let start = self.offset + self.curpos;
            let end = self.offset + self.curpos + self.k;
            self.curpos += self.k;
            Some((
                self.seq[start..end].to_vec(),
                (start, end-1)
            ))
        }
    }
}

// Classifies Kmers as gff3 entries...
#[derive(Debug, PartialEq)]
pub struct Gff3Kmers {
    pub kmers: Vec<Vec<u8>>,
    pub id: String,
    pub classification: Vec<Vec<u8>>,
}

pub struct Gff3KmersGenerator {
    kmercoords_window_generator: KmerCoordsWindowGenerator,
    intervals: intervals::IntervalMap<Vec<u8>>,
    pub types: Vec<String>,

    // TODO: Put alphabet here so we can train proteins as well...
}

impl Gff3KmersGenerator {
    pub fn new(
        gff3file: String,
        generator: KmerCoordsWindowGenerator) 
    -> Gff3KmersGenerator {
        let (intervals, types) = gff3::get_gff3_intervals(gff3file);

        Gff3KmersGenerator {
            kmercoords_window_generator: generator,
            intervals,
            types,
        }
    }
}

impl Iterator for Gff3KmersGenerator {
    type Item = Gff3Kmers;

    fn next(&mut self) -> Option<Gff3Kmers> {
        let next_item = match self.kmercoords_window_generator.next() {
            Some(x) => x,
            None    => return None
        };

        let KmerCoordsWindow { kmers, id, coords } = next_item;

        let mut classification = Vec::new();

        let blank: Vec<u8> = vec![0; self.types.len()];

        for (n, _) in kmers.iter().enumerate() {
            let found = 
              self.intervals
                .landmarks
                .get(&id)
                .expect("GFF3 File Error")
                .find(coords[n].0 as u32, coords[n].1 as u32);

            let mut kmer_classification = blank.clone();

            for x in found {
                for (i, m) in (*x).val.iter().enumerate() {
                    kmer_classification[i] += m;
                }
            }

            classification.push(kmer_classification);
        }

        Some(Gff3Kmers { kmers, id, classification })
    }
}

/*pub struct DiscriminatorMaskedA2T {
    pub kmers: Vec<Vec<u8>>,
    pub id: String,
    pub taxons: Vec<usize>,
    pub taxon: usize,
    pub truth: Vec<u8>,
} */


#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::prelude::*;
    use super::*;

    #[test]
    pub fn test_kmers_iter() {
        let mut kmers = Kmers::new(b"ACTGACTGACTGACTG".to_vec(), 3, 0);

        let k1 = kmers.next().expect("Unable to get Kmer");
        assert!("ACT".to_string() == std::str::from_utf8(&k1.0).expect("Unable to convert from Vec<u8>"));
        
        let k2 = kmers.next().expect("Unable to get Kmer");
        assert!("GAC".to_string() == std::str::from_utf8(&k2.0).expect("Unable to convert from Vec<u8>"));
        assert!((3,5) == k2.1);
    }


    #[test]
    pub fn test_kmer_window_generator() {
        let mut kmers = KmerWindowGenerator::new("test_data/test.sfasta".to_string(), 3, 3, 0, false);
        let first = kmers.next().expect("Unable to get KmerWindow");

        assert!(first.kmers[0] == b"ACT");

        let mut kmers = KmerWindowGenerator::new("test_data/test.sfasta".to_string(), 3, 3, 0, false);

        let mut skipped = kmers.skip(2).next().expect("Unable to skip ahead");
        assert!(skipped.kmers[0] == b"NAC");
    }

    #[test]
    pub fn test_kmer_coords_window_generator() {
        let mut kmers = KmerCoordsWindowGenerator::new("test_data/test.sfasta".to_string(), 3, 3, 0, false);
        let first = kmers.next().expect("Unable to get KmerWindow");
        println!("{:#?}", first.coords);

        assert!(first.coords == [(0,2),(3,5),(6,8)]);
        assert!(first.kmers[0] == b"ACT");

        let mut kmers = KmerCoordsWindowGenerator::new("test_data/test.sfasta".to_string(), 3, 3, 0, false);

        let mut skipped = kmers.skip(2).next().expect("Unable to skip ahead");
        assert!(skipped.kmers[0] == b"NAC");
        println!("{:#?}", skipped.coords);
        assert!(skipped.coords == [(38,40),(41,43),(44,46)]);
    }
}