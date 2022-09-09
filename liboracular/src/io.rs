pub use libsfasta as sfasta;
pub use sfasta::prelude::*;

pub struct SequenceSplitter3N {
    sequences: Box<dyn Iterator<Item = Sequence> + Send>,
    curseq: Sequence,
    curpos: usize,
    curlen: usize,
}

impl SequenceSplitter3N {
    pub fn new(mut sequences: Box<dyn Iterator<Item = Sequence> + Send>) -> SequenceSplitter3N {
        let curseq = match sequences.next() {
            Some(x) => x,
            None => panic!("File is empty!"),
        };

        // println!("{:#?}", std::str::from_utf8(&curseq.seq).unwrap());

        let curlen = curseq.len();

        SequenceSplitter3N {
            sequences,
            curseq,
            curpos: 0,
            curlen,
        }
    }

    fn next_seq(&mut self) -> bool {
        self.curpos = 0;
        let curseq = match self.sequences.next() {
            Some(x) => x,
            None => {
                // println!("Stop1");
                return false;
            }
        };
        self.curlen = curseq.len();
        self.curseq = curseq;
        return true;
    }
}

impl Iterator for SequenceSplitter3N {
    type Item = Sequence;

    fn next(&mut self) -> Option<Sequence> {
        if self.curlen == self.curpos && !self.next_seq() {
            return None;
        }

        let startloc;
        let mut endloc;

        loop {
            if bytecount::count(&self.curseq.sequence.as_ref().unwrap(), b'N') < 3 {
                startloc = 0;
                endloc = self.curlen;
            } else {
                startloc = self.curpos
                    + match self.curseq.sequence.as_ref().unwrap()[self.curpos..]
                        .windows(3)
                        .enumerate()
                        .filter(|(_y, x)| bytecount::count(&x, b'N') < 3)
                        .map(|(y, _x)| y)
                        .nth(0)
                    {
                        Some(x) => x,
                        None => {
                            if !self.next_seq() {
                                return None;
                            }
                            continue;
                        }
                    };

                endloc = startloc
                    + 1
                    + self.curseq.sequence.as_ref().unwrap()[startloc + 1..]
                        .windows(3)
                        .enumerate()
                        .filter(|(_y, x)| bytecount::count(&x, b'N') == 3)
                        .map(|(y, _x)| y)
                        .nth(0)
                        .unwrap_or(self.curlen - startloc);
            }

            if endloc > self.curseq.len() {
                endloc = self.curseq.len();
            }

            break;
        }

        // println!("{} {} {}", startloc, endloc, self.curseq.seq.len());

        self.curpos = endloc;

        Some(Sequence {
            id: self.curseq.id.clone(),
            sequence: Some(self.curseq.sequence.as_ref().unwrap()[startloc..endloc].to_vec()),
            header: None,
            scores: None,
            offset: startloc,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    pub fn test_iosequences() {
        convert_fasta_file("test_data/test.fna", "test_data/test_sequences.sfasta");

        let mut sequences = Box::new(Sequences::from_file("test_data/test_sequences.sfasta"));
        let sequence = sequences.next().unwrap();
        println!("{:#?}", sequence);
        assert!(sequence.id.as_ref().unwrap() == "test");
        assert!(sequence.len() == 670);
    }

    #[test]
    #[should_panic]
    pub fn test_empty() {
        let sequences = Box::new(Sequences::from_file("test_data/empty.sfasta"));
        Box::new(SequenceSplitter3N::new(sequences));
    }

    #[test]
    pub fn test_3n_splitter() {
        convert_fasta_file("test_data/test.fna", "test_data/test_3n_splitter.sfasta");
        let sequences = Box::new(Sequences::from_file("test_data/test_3n_splitter.sfasta"));
        let count = sequences.count();
        assert!(count == 1);

        let sequences = Box::new(Sequences::from_file("test_data/test_3n_splitter.sfasta"));
        let mut sequences = Box::new(SequenceSplitter3N::new(sequences));
        let x = sequences.next().unwrap();
        println!("{}", x.len());
        println!(
            "{}",
            std::str::from_utf8(&x.sequence.as_ref().unwrap()).unwrap()
        );
        assert!(x.len() == 20);

        let x = sequences.next().unwrap();
        println!("{}", x.len());
        println!(
            "{}",
            std::str::from_utf8(&x.sequence.as_ref().unwrap()).unwrap()
        );
        assert!(x.len() == 50);

        let x = sequences.next().unwrap();
        println!("{}", x.len());
        println!(
            "{}",
            std::str::from_utf8(&x.sequence.as_ref().unwrap()).unwrap()
        );
        assert!(x.len() == 50);

        let x = sequences.next().unwrap();
        println!("{}", x.len());
        println!(
            "{}",
            std::str::from_utf8(&x.sequence.as_ref().unwrap()).unwrap()
        );
        assert!(x.len() == 50);

        let x = sequences.next().unwrap();
        println!("{}", x.len());
        println!(
            "{}",
            std::str::from_utf8(&x.sequence.as_ref().unwrap()).unwrap()
        );
        assert!(x.len() == 50);

        let mut therest: Vec<Sequence> = sequences.collect();
        let x = therest.pop().unwrap();
        println!("{}", x.len());
        println!(
            "{}",
            std::str::from_utf8(&x.sequence.as_ref().unwrap()).unwrap()
        );
        assert!(x.len() == 30);

        convert_fasta_file(
            "test_data/test_multiple.fna",
            "test_data/test_3n_splitter.sfasta",
        );

        let mut sequences = Box::new(Sequences::from_file("test_data/test_3n_splitter.sfasta"));
        sequences.set_mode(SeqMode::Random);
        let mut sequences = Box::new(SequenceSplitter3N::new(sequences));

        let count = sequences.count();
        println!("Sequences after 3N Split: {}", count);
        assert!(count == 216);

        let sequences = Box::new(Sequences::from_file("test_data/test_3n_splitter.sfasta"));
        let mut sequences = Box::new(SequenceSplitter3N::new(sequences));
        sequences.next();
        // println!("Coords: {:#?}", sequences.coords[0]);
        // assert!(sequences.coords[0] == (38, 86));
    }

    #[test]
    pub fn test_sequences_impl() {
        convert_fasta_file(
            "test_data/test_multiple.fna",
            "test_data/test_sequences_impl.sfasta",
        );
        let seqs = Box::new(Sequences::from_file("test_data/test_sequences_impl.sfasta"));
        let count = seqs.count();
        println!("Sequences Impl Count {}", count);
        assert!(count == 8);
    }

    /*
    TODO: Fix this, but upstream in libsfasta
    #[test]
    #[should_panic(expected = "Couldn't open test_data/empty.fna.sfasta")]
    pub fn test_3n_splitter_empty() {
        convert_fasta_file("test_data/empty.fna", "test_data/empty.fna.sfasta");
        let sequences = Box::new(Sequences::from_file("test_data/empty.fna.sfasta"));
        let _ = Box::new(SequenceSplitter3N::new(sequences));
    }
    */
}
