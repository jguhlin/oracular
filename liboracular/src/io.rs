use std::fs::File;
use std::io::{BufReader, Read};

use serde::{Deserialize, Serialize};

use std::collections::VecDeque;

use crate::io;

pub struct Sequences {
    reader: Box<dyn Read + Send>,
}

#[derive(PartialEq, Serialize, Deserialize, Clone, Debug)]
pub struct Sequence {
    pub seq: Vec<u8>,
    pub id: String,
    pub location: usize,
    pub end: usize,
}

pub struct SequenceSplitter3N {
    sequences: Box<dyn Iterator<Item = io::Sequence> + Send>,
    curseq: Sequence,
    coords: VecDeque<(usize, usize)>,
}

impl SequenceSplitter3N {
    pub fn new(mut sequences: Box<dyn Iterator<Item = io::Sequence> + Send>) -> SequenceSplitter3N {
        let curseq = match sequences.next() {
            Some(x) => x,
            None => panic!("File is empty!"),
        };

        let coords: VecDeque<(usize, usize)>;
        coords = crate::utils::get_good_sequence_coords(&curseq.seq)
            .into_iter()
            .collect();

        SequenceSplitter3N {
            sequences,
            curseq,
            coords,
        }
    }
}

// TODO: This is the right place to do this, but I feel it's happening somewhere
// else and wasting CPU cycles...
impl Sequence {
    pub fn make_uppercase(&mut self) {
        self.seq.make_ascii_uppercase();
    }
}

impl Iterator for SequenceSplitter3N {
    type Item = Sequence;

    #[inline]
    fn next(&mut self) -> Option<Sequence> {
        let coords = match self.coords.pop_front() {
            Some(x) => x,
            None => {
                let mut x = None;
                while x == None {
                    let curseq = match self.sequences.next() {
                        Some(x) => x,
                        None => return None,
                    };

                    let coords: VecDeque<(usize, usize)>;
                    coords = crate::utils::get_good_sequence_coords(&curseq.seq)
                        .into_iter()
                        .collect();

                    self.curseq = curseq;
                    self.coords = coords;
                    x = match self.coords.pop_front() {
                        Some(x) => Some(x),
                        None => None,
                    }
                }

                x.unwrap()
            }
        };

        Some(Sequence {
            id: self.curseq.id.clone(),
            seq: self.curseq.seq[coords.0..coords.1].to_vec(),
            location: coords.0 as usize,
            end: coords.1 as usize,
        })
    }
}

fn open_file(filename: String) -> Box<dyn Read + Send> {
    let file = match File::open(&filename) {
        Err(why) => panic!("Couldn't open {}: {}", filename, why.to_string()),
        Ok(file) => file,
    };

    let file = BufReader::with_capacity(64 * 1024 * 1024, file);

    let fasta: Box<dyn Read + Send> = if filename.ends_with("gz") {
        Box::new(flate2::read::GzDecoder::new(file))
    } else if filename.ends_with("snappy") || filename.ends_with("sz") {
        Box::new(snap::read::FrameDecoder::new(file))
    } else {
        Box::new(file)
    };

    let reader = BufReader::with_capacity(32 * 1024 * 1024, fasta);
    Box::new(reader)
}

impl Iterator for Sequences {
    type Item = Sequence;

    fn next(&mut self) -> Option<Sequence> {
        let seq: Sequence = match bincode::deserialize_from(&mut self.reader) {
            Ok(x) => x,
            Err(_) => return None,
        };
        Some(seq)
    }
}

impl Sequences {
    pub fn new(filename: String) -> Sequences {
        Sequences {
            reader: open_file(filename),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sfasta;

    #[test]
    pub fn test_sequences() {
        let mut sequences = Box::new(sfasta::Sequences::new("test_data/test.sfasta"));
        let sequence = sequences.next().unwrap();
        assert!(sequence.id == "test");
        assert!(sequence.end == 669);
        assert!(sequence.location == 0);
        assert!(sequence.seq.len() == 669);
    }
    /*
    #[test]
    pub fn test_regular_fasta_file() {
        let mut sequences = Sequences::new("test_data/test.fna".to_string());
        let sequence = sequences.next().unwrap();
        assert!(sequence.id == "test");
        assert!(sequence.end == 669);
        assert!(sequence.location == 0);
        assert!(sequence.seq.len() == 669);




    }*/

    #[test]
    #[should_panic]
    pub fn test_empty() {
        let sequences = Box::new(sfasta::Sequences::new("test_data/empty.sfasta"));
        Box::new(io::SequenceSplitter3N::new(sequences));
    }

    #[test]
    pub fn test_3n_splitter() {
        let sequences = Box::new(sfasta::Sequences::new("test_data/test.sfasta"));
        let sequences = Box::new(io::SequenceSplitter3N::new(sequences));
        assert!(sequences.coords[0] == (0, 19));

        sfasta::convert_fasta_file(
            "test_data/test_multiple.fna",
            "test_data/test_multiple.sfasta",
        );
        let sequences = Box::new(sfasta::Sequences::new("test_data/test_multiple.sfasta"));
        let mut sequences = Box::new(io::SequenceSplitter3N::new(sequences));
        sequences.next();
        println!("Coords: {:#?}", sequences.coords[0]);
        assert!(sequences.coords[0] == (38, 79));

        sequences.next().unwrap();
        sequences.next().unwrap();
        sequences.next().unwrap();
        sequences.next().unwrap();
        let x = sequences.next().unwrap();

        // println!("{:#?}", x);
        assert!(x.id == "test2");
        assert!(x.seq[0..5] == [78, 65, 67, 84, 71]);
        assert!(x.location == 105);
        assert!(x.end == 153);
    }

    #[test]
    pub fn test_sequences_impl() {
        Sequences::new("test_data/test.fna".to_string());
    }

    #[test]
    #[should_panic(expected = "Couldn't open test_data/empty.fna.sfasta")]
    pub fn test_3n_splitter_empty() {
        let sequences = Box::new(sfasta::Sequences::new("test_data/empty.fna"));
        let _ = Box::new(io::SequenceSplitter3N::new(sequences));
    }
}
