use rand::prelude::*;
use rand::Rng;
use std::fmt;

use crate::gff3;
use crate::intervals;
use crate::io;
use crate::sfasta;
use crate::utils;

type Kmer = Vec<u8>;
type Coords = (usize, usize);

pub struct KmerCoords {
    pub kmer: Vec<u8>,
    pub coords: Coords,
}

impl fmt::Display for KmerCoords {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}: {}-{}",
            std::str::from_utf8(&self.kmer).expect("Unable to convert Kmer to Str"),
            self.coords.0,
            self.coords.1,
        )
    }
}

impl fmt::Debug for KmerCoords {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}: {}-{}",
            std::str::from_utf8(&self.kmer).expect("Unable to convert Kmer to Str"),
            self.coords.0,
            self.coords.1,
        )
    }
}

impl std::cmp::PartialEq<&[u8]> for KmerCoords {
    fn eq(&self, other: &&[u8]) -> bool {
        self.kmer == *other
    }
}

#[derive(Debug, Clone)]
pub struct KmerWindow {
    pub kmers: Vec<Kmer>,
    pub id: String,
    pub rc: bool,
}

pub struct KmerWindowGenerator {
    sequences: Box<dyn Iterator<Item = io::Sequence> + Send>,
    window_size: usize,
    k: usize,
    kmer_generator: Kmers,
    needed_sequence: usize,
    curseq: io::Sequence,
    offset: usize,
    rc: bool,
}

#[derive(Debug, PartialEq)]
pub struct DiscriminatorMasked {
    pub kmers: Vec<Vec<u8>>,
    pub id: String,
    pub truth: Vec<u8>,
}

pub struct DiscriminatorMaskedGenerator {
    kmer_window_generator: KmerWindowGenerator,
    replacement_pct: f32,
    k: usize,
    // TODO: Put alphabet here so we can train proteins as well...
}

impl DiscriminatorMaskedGenerator {
    pub fn new(
        replacement_pct: f32,
        k: usize,
        generator: KmerWindowGenerator,
    ) -> DiscriminatorMaskedGenerator {
        // let rng = Arc::new(RwLock::new(rand::thread_rng()));

        DiscriminatorMaskedGenerator {
            kmer_window_generator: generator,
            replacement_pct,
            k,
        }
    }
}

impl Iterator for DiscriminatorMaskedGenerator {
    type Item = DiscriminatorMasked;

    fn next(&mut self) -> Option<DiscriminatorMasked> {
        let next_item = match self.kmer_window_generator.next() {
            Some(x) => x,
            None => return None,
        };

        let KmerWindow {
            mut kmers,
            id,
            rc: _,
        } = next_item;

        // TODO: Make switchable, so we can train protein sequences
        // ~2% chance of an N
        let choices = [(b'A', 48), (b'C', 48), (b'T', 48), (b'G', 48), (b'N', 4)];

        let mut truth: Vec<u8> = Vec::with_capacity(kmers.len());

        let mut rng = rand::thread_rng();

        for kmer in kmers.iter_mut() {
            if rng.gen::<f32>() < self.replacement_pct {
                let mut new_kmer: Vec<u8> = Vec::with_capacity(self.k);
                for _i in 0..self.k {
                    new_kmer.push(choices.choose_weighted(&mut rng, |item| item.1).unwrap().0);
                }
                *kmer = new_kmer;
                truth.push(0);
            } else {
                truth.push(1);
            }
        }

        // return Some(DiscriminatorMasked { kmers, id, taxons, taxon, truth })
        Some(DiscriminatorMasked { kmers, id, truth })
    }
}

impl KmerWindowGenerator {
    pub fn new(
        filename: &str,
        k: usize,
        window_size: usize,
        offset: usize,
        rc: bool,
        rand: bool, // Mostly for debugging, but also for synteny plots and such...
    ) -> KmerWindowGenerator {
        let mut sequences = Box::new(sfasta::Sequences::new(filename));
        if rand {
            sequences.set_mode(sfasta::SeqMode::Random);
        }

        let mut sequences = Box::new(io::SequenceSplitter3N::new(sequences));
        let mut curseq = match sequences.next() {
            Some(x) => x,
            None => panic!("File is empty or invalid format!"),
        };

        if rc {
            let io::Sequence {
                id,
                mut seq,
                location,
                end,
            } = curseq;

            utils::complement_nucleotides(&mut seq);
            seq.reverse();

            // TODO: How to deal with rc?

            curseq = io::Sequence {
                id,
                seq,
                location,
                end,
            };
        }

        let kmer_generator = Kmers::new(curseq.seq.clone(), k, offset, rc);
        // TODO: Should be able to handle between min_seq and max_seq
        // So instead of only 20 or 30 kmers at a time, get 10 - 100 (prefer more)
        //let needed_sequence = k * window_size;
        let needed_sequence = k;
        // HERE! Need to return small sequences unless < window_size, and let the generators deal with removing them or not...

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
        // While instead of if, because if we get a too short sequence we should skip
        // it...
        // and HERE! Need to return small sequences unless < window_size, and let the generators deal with removing them or not...
        while (self.kmer_generator.len - self.kmer_generator.curpos) <= self.needed_sequence {
            let mut curseq: io::Sequence = match self.sequences.next() {
                Some(x) => x,
                None => return None, // That's it... no more!
            };

            if self.rc {
                let io::Sequence {
                    id,
                    mut seq,
                    location,
                    end,
                } = curseq;

                utils::complement_nucleotides(&mut seq);
                seq.reverse();

                curseq = io::Sequence {
                    id,
                    seq,
                    location,
                    end,
                };
            }

            self.curseq = curseq.clone();
            let kmer_generator = Kmers::new(curseq.seq.clone(), self.k, self.offset, self.rc);
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
                        panic!("Invalid Kmer");
                    }
                }
                None => return None,
            };
        }

        Some(KmerWindow {
            kmers,
            id: self.curseq.id.clone(),
            rc: self.rc,
        })
    }
}

#[derive(Debug)]
pub struct KmerCoordsWindow {
    pub kmers: Vec<Kmer>,
    pub coords: Vec<Coords>,
    pub id: String,
    pub rc: bool,
}

pub struct KmerCoordsWindowIter {
    sequences: Box<dyn Iterator<Item = io::Sequence> + Send>,
    window_size: usize,
    k: usize,
    kmer_generator: Kmers,
    needed_sequence: usize,
    curseq: io::Sequence,
    offset: usize,
    rc: bool,
}

impl KmerCoordsWindowIter {
    pub fn new(
        filename: &str,
        k: usize,
        window_size: usize,
        offset: usize,
        rc: bool,
        rand: bool,
    ) -> KmerCoordsWindowIter {
        let mut sequences = Box::new(sfasta::Sequences::new(filename));
        if rand {
            sequences.set_mode(sfasta::SeqMode::Random);
        }

        let mut sequences = Box::new(io::SequenceSplitter3N::new(sequences));
        let mut curseq = match sequences.next() {
            Some(x) => x,
            None => panic!("File is empty or invalid format!"),
        };

        if rc {
            let io::Sequence {
                id,
                mut seq,
                location,
                end,
            } = curseq;

            utils::complement_nucleotides(&mut seq);
            seq.reverse();

            curseq = io::Sequence {
                id,
                seq,
                location,
                end,
            };
        }

        let kmer_generator = Kmers::new(curseq.seq.clone(), k, offset, rc);
        let needed_sequence = k * window_size;

        KmerCoordsWindowIter {
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

impl Iterator for KmerCoordsWindowIter {
    type Item = KmerCoordsWindow;

    fn next(&mut self) -> Option<KmerCoordsWindow> {
        // While instead of if, because if we get a too short sequence we should skip
        // it...
        while (self.kmer_generator.len - self.kmer_generator.curpos) <= self.needed_sequence {
            let mut curseq: io::Sequence = match self.sequences.next() {
                Some(x) => x,
                None => return None, // That's it... no more!
            };

            if self.rc {
                let io::Sequence {
                    id,
                    mut seq,
                    location,
                    end,
                } = curseq;

                utils::complement_nucleotides(&mut seq);
                seq.reverse();

                curseq = io::Sequence {
                    id,
                    seq,
                    location,
                    end,
                };
            }

            self.curseq = curseq.clone();
            let kmer_generator = Kmers::new(curseq.seq.clone(), self.k, self.offset, self.rc);
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
                }
                None => return None,
            };
        }

        /*        if self.rc {
            coords = coords.iter().map(
                |(x, y)|
                (self.curseq.location+self.curseq.end-y,
                 self.curseq.location+self.curseq.end-x-1,)).collect();
        } else { */
        coords = coords
            .iter()
            .map(|(x, y)| (x + self.curseq.location, y + self.curseq.location))
            .collect();
        //}

        Some(KmerCoordsWindow {
            kmers,
            id: self.curseq.id.clone(),
            coords,
            rc: self.rc,
        })
    }
}

pub struct Kmers {
    seq: Vec<u8>,
    k: usize,
    offset: usize,
    curpos: usize,
    len: usize,
    rc: bool,
}

impl Kmers {
    pub fn new(seq: Vec<u8>, k: usize, offset: usize, rc: bool) -> Kmers {
        let len = seq.len();

        Kmers {
            seq,
            k,
            offset,
            len,
            curpos: 0,
            rc,
        }
    }
}

impl Iterator for Kmers {
    type Item = (Vec<u8>, Coords);

    fn next(&mut self) -> Option<(Vec<u8>, Coords)> {
        if (self.offset + self.curpos + self.k) >= self.len {
            None
        } else {
            let start = self.offset + self.curpos;
            let end = self.offset + self.curpos + self.k;
            self.curpos += self.k;
            let coords;
            if self.rc {
                // Start counting from the BACK of the sequence
                // Sequence already represents RC, but coords do not...
                coords = (self.len - end, self.len - start - 1)
            } else {
                coords = (start, end - 1)
            }

            Some((self.seq[start..end].to_vec(), coords))
        }
    }
}

// Classifies Kmers as gff3 entries...
#[derive(Debug, PartialEq)]
pub struct Gff3Kmers {
    pub kmers: Vec<Vec<u8>>,
    pub id: String,
    pub classifications: Vec<Vec<u8>>,
    pub coords: Vec<Coords>,
    pub rc: bool,
}

pub struct Gff3KmersIter {
    kmercoords_window_iter: KmerCoordsWindowIter,
    intervals: intervals::IntervalMap<Vec<u8>>,
    pub types: Vec<String>,
    pub k: usize,
}

impl Gff3KmersIter {
    pub fn new(gff3file: &str, generator: KmerCoordsWindowIter, k: usize) -> Gff3KmersIter {
        let (intervals, types) = gff3::get_gff3_intervals(gff3file);

        Gff3KmersIter {
            kmercoords_window_iter: generator,
            intervals,
            types,
            k,
        }
    }
}

impl Iterator for Gff3KmersIter {
    type Item = Gff3Kmers;

    fn next(&mut self) -> Option<Gff3Kmers> {
        let next_item = match self.kmercoords_window_iter.next() {
            Some(x) => x,
            None => return None,
        };

        let KmerCoordsWindow {
            kmers,
            id,
            coords,
            rc,
        } = next_item;

        let mut classifications = Vec::new();

        let blank: Vec<u8> = vec![0; self.types.len()];

        for (n, _) in kmers.iter().enumerate() {
            let found = match self.intervals.landmarks.get(&id) {
                Some(x) => Some(x.find(coords[n].0 as u32, coords[n].1 as u32)),
                None => None,
            };

            let mut kmer_classification = blank.clone();

            if let Some(found) = found {
                for x in found {
                    let lower: u32 = std::cmp::max(x.start as u32, coords[n].0 as u32);
                    let upper: u32 = std::cmp::min(x.stop as u32, coords[n].1 as u32);

                    let len = upper - lower;
                    if len as f32 > (self.k as f32 * 0.8) {
                        // TODO: Make overlap (0.8) configurable...
                        // Intervals are a set of [0 0 0 0 1 1 1 1 1] etc... for when classes are
                        // false or true...
                        for (i, m) in (*x).val.iter().enumerate() {
                            if *m == 1 {
                                kmer_classification[i] = 1;
                            }
                        }
                    }
                }
            }

            classifications.push(kmer_classification);
        }

        Some(Gff3Kmers {
            kmers,
            id,
            classifications,
            coords,
            rc,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    pub fn test_kmers_iter() {
        let mut kmers = Kmers::new(b"ACTGACTGACTGACTG".to_vec(), 3, 0, false);

        let k1 = kmers.next().expect("Unable to get Kmer");
        assert!("ACT" == std::str::from_utf8(&k1.0).expect("Unable to convert from Vec<u8>"));

        let k2 = kmers.next().expect("Unable to get Kmer");
        assert!("GAC" == std::str::from_utf8(&k2.0).expect("Unable to convert from Vec<u8>"));
        assert!((3, 5) == k2.1);
    }

    #[test]
    pub fn test_kmer_window_generator() {
        crate::sfasta::convert_fasta_file(
            "test_data/test.fna",
            "test_data/test_kmer_coords_window_generator.sfasta",
        );
        let mut kmers = KmerWindowGenerator::new(
            "test_data/test_kmer_coords_window_generator.sfasta",
            3,
            3,
            0,
            false,
            false,
        );
        let first = kmers.next().expect("Unable to get KmerWindow");

        assert!(first.kmers[0] == b"ACT");

        let mut kmers = KmerWindowGenerator::new("test_data/test.sfasta", 3, 3, 0, false, false);

        let skipped = kmers.nth(2).expect("Unable to skip ahead");
        assert!(skipped.kmers[0] == b"NAC");
    }

    #[test]
    pub fn test_kmers_rc_coords() {
        let mut kmers = Kmers::new(b"ACTGACTGACTGACTG".to_vec(), 3, 0, true);

        let k1 = kmers.next().expect("Unable to get Kmer");
        assert!(k1.1 == (13, 15));
        assert!("ACT" == std::str::from_utf8(&k1.0).expect("Unable to convert from Vec<u8>"));

        let kmers = Kmers::new(b"ACTGACTGACTGACTG".to_vec(), 3, 0, true);
        let coords: Vec<_> = kmers.map(|x| x.1).collect();

        println!("{:#?}", coords);
    }

    #[test]
    pub fn test_kmer_coords_window_generator() {
        let mut kmers = KmerCoordsWindowIter::new("test_data/test.sfasta", 3, 3, 0, false, false);
        let first = kmers.next().expect("Unable to get KmerWindow");

        assert!(first.coords == [(0, 2), (3, 5), (6, 8)]);
        assert!(first.kmers[0] == b"ACT");

        let mut kmers = KmerCoordsWindowIter::new("test_data/test.sfasta", 3, 3, 0, false, false);

        let skipped = kmers.nth(2).expect("Unable to skip ahead");
        // println!("{:#?}", skipped.kmers[0]);
        assert!(skipped.kmers[0] == b"NAC");
        // println!("{:#?}", skipped.coords);
        assert!(skipped.coords == [(38, 40), (41, 43), (44, 46)]);

        // Have to get the right coords for RC
        let kmers = KmerCoordsWindowIter::new("test_data/test.sfasta", 8, 2, 0, true, false);

        let coords: Vec<_> = kmers.map(|x| x.coords).collect();
        println!("{:#?}", coords[0]);
        println!("{:#?}", coords[1]);
        assert!(coords[0] == [(11, 18), (3, 10)]);
        assert!(coords[1] == [(78, 85), (70, 77)]);
    }
}
