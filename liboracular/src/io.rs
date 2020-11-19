use std::fs::File;
use std::io::{BufReader, Read};

use serde::{Deserialize, Serialize};

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
    curpos: usize,
    curlen: usize,
}

impl SequenceSplitter3N {
    pub fn new(mut sequences: Box<dyn Iterator<Item = io::Sequence> + Send>) -> SequenceSplitter3N {
        let curseq = match sequences.next() {
            Some(x) => x,
            None => panic!("File is empty!"),
        };

        // println!("{:#?}", std::str::from_utf8(&curseq.seq).unwrap());

        let curlen = curseq.seq.len();

        SequenceSplitter3N {
            sequences,
            curseq,
            curpos: 0,
            curlen,
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

    fn next(&mut self) -> Option<Sequence> {
        if self.curlen == self.curpos {
            self.curpos = 0;
            let curseq = match self.sequences.next() {
                Some(x) => x,
                None => {
                    // println!("Stop1");
                    return None;
                }
            };
            self.curlen = curseq.seq.len();
            self.curseq = curseq;
        }

        let startloc;
        let mut endloc;

        if bytecount::count(&self.curseq.seq, b'N') < 3 {
            startloc = 0;
            endloc = self.curlen;
        } else {
            startloc = self.curpos
                + match self.curseq.seq[self.curpos..]
                    .windows(3)
                    .enumerate()
                    .filter(|(_y, x)| bytecount::count(&x, b'N') < 3)
                    .map(|(y, _x)| y)
                    .nth(0)
                {
                    Some(x) => x,
                    None => {
                        // println!("Stop2");
                        return None;
                    }
                };

            endloc = startloc
                + 1
                + self.curseq.seq[startloc + 1..]
                    .windows(3)
                    .enumerate()
                    .filter(|(_y, x)| bytecount::count(&x, b'N') == 3)
                    .map(|(y, _x)| y)
                    .nth(0)
                    .unwrap_or(self.curlen - startloc);
        }

        if endloc > self.curseq.seq.len() {
            endloc = self.curseq.seq.len();
        }

        // println!("{} {} {}", startloc, endloc, self.curseq.seq.len());

        self.curpos = endloc;

        Some(Sequence {
            id: self.curseq.id.clone(),
            seq: self.curseq.seq[startloc..endloc].to_vec(),
            location: startloc as usize,
            end: endloc as usize,
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
            Err(_) => {
                println!("SeqStop");
                return None;
            }
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
    pub fn test_iosequences() {
        sfasta::clear_idxcache();
        sfasta::convert_fasta_file("test_data/test.fna", "test_data/test_sequences.sfasta");
        let mut sequences = Box::new(sfasta::Sequences::new("test_data/test_sequences.sfasta"));
        let sequence = sequences.next().unwrap();
        println!("{:#?}", sequence);
        assert!(sequence.id == "test");
        assert!(sequence.end == 670);
        assert!(sequence.location == 0);
        assert!(sequence.seq.len() == 670);
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
        sfasta::clear_idxcache();
        sfasta::convert_fasta_file("test_data/test.fna", "test_data/test_3n_splitter.sfasta");
        let sequences = Box::new(sfasta::Sequences::new("test_data/test_3n_splitter.sfasta"));
        let count = sequences.count();
        assert!(count == 1);

        let sequences = Box::new(sfasta::Sequences::new("test_data/test_3n_splitter.sfasta"));
        let mut sequences = Box::new(io::SequenceSplitter3N::new(sequences));
        let x = sequences.next().unwrap();
        println!("{}", x.seq.len());
        println!("{}", std::str::from_utf8(&x.seq).unwrap());
        assert!(x.seq.len() == 20);

        let x = sequences.next().unwrap();
        println!("{}", x.seq.len());
        println!("{}", std::str::from_utf8(&x.seq).unwrap());
        assert!(x.seq.len() == 50);

        let x = sequences.next().unwrap();
        println!("{}", x.seq.len());
        println!("{}", std::str::from_utf8(&x.seq).unwrap());
        assert!(x.seq.len() == 50);

        let x = sequences.next().unwrap();
        println!("{}", x.seq.len());
        println!("{}", std::str::from_utf8(&x.seq).unwrap());
        assert!(x.seq.len() == 50);

        let x = sequences.next().unwrap();
        println!("{}", x.seq.len());
        println!("{}", std::str::from_utf8(&x.seq).unwrap());
        assert!(x.seq.len() == 50);

        let mut therest: Vec<Sequence> = sequences.collect();
        let x = therest.pop().unwrap();
        println!("{}", x.seq.len());
        println!("{}", std::str::from_utf8(&x.seq).unwrap());
        assert!(x.seq.len() == 30);

        sfasta::clear_idxcache();
        sfasta::convert_fasta_file(
            "test_data/test_multiple.fna",
            "test_data/test_3n_splitter.sfasta",
        );

        let mut sequences = Box::new(sfasta::Sequences::new("test_data/test_3n_splitter.sfasta"));
        sequences.set_mode(sfasta::SeqMode::Random);
        let mut sequences = Box::new(io::SequenceSplitter3N::new(sequences));

        let count = sequences.count();
        println!("Sequences after 3N Split: {}", count);
        assert!(count == 216);

        let sequences = Box::new(sfasta::Sequences::new("test_data/test_3n_splitter.sfasta"));
        let mut sequences = Box::new(io::SequenceSplitter3N::new(sequences));
        sequences.next();
        // println!("Coords: {:#?}", sequences.coords[0]);
        // assert!(sequences.coords[0] == (38, 86));
    }

    #[test]
    pub fn test_sequences_impl() {
        sfasta::clear_idxcache();
        sfasta::convert_fasta_file(
            "test_data/test_multiple.fna",
            "test_data/test_sequences_impl.sfasta",
        );
        let seqs = Box::new(sfasta::Sequences::new(
            "test_data/test_sequences_impl.sfasta",
        ));
        let count = seqs.count();
        println!("Sequences Impl Count {}", count);
        assert!(count == 8);
    }

    #[test]
    #[should_panic(expected = "Couldn't open test_data/empty.fna.sfasta")]
    pub fn test_3n_splitter_empty() {
        let sequences = Box::new(sfasta::Sequences::new("test_data/empty.fna"));
        let _ = Box::new(io::SequenceSplitter3N::new(sequences));
    }
}
