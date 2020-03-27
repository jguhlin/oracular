use std::fs::File;
use std::io::{BufReader, Read, BufRead};

use indicatif::ProgressBar;
use indicatif::ProgressStyle;

use serde::{Serialize, Deserialize};

#[derive(PartialEq, Serialize, Deserialize)]
pub struct Sequence {
    pub seq: Vec<u8>,
    pub id:  String,
}

pub struct Sequences {
    reader: Box<dyn Read>
}

fn open_file_with_progress_bar(filename: String) -> (Box<dyn Read>, ProgressBar)
{
    let file = match File::open(&filename) {
        Err(why) => panic!("Couldn't open {}: {}", filename, why.to_string()),
        Ok(file) => file,
    };

    let pb = ProgressBar::new(file.metadata().unwrap().len());
    pb.set_style(ProgressStyle::default_bar()
                    .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>5}/{len:5} {eta_precise} {msg}")
                    .progress_chars("█▇▆▅▄▃▂▁  "));

    let file = BufReader::with_capacity(64 * 1024 * 1024, pb.wrap_read(file));

    let fasta: Box<dyn Read> = if filename.ends_with("gz") {
        Box::new(flate2::read::GzDecoder::new(file))
    } else if filename.ends_with("snappy") || filename.ends_with("sz") {
        Box::new(snap::read::FrameDecoder::new(file))
    } else {
        Box::new(file)
    };

    let reader = BufReader::with_capacity(32 * 1024 * 1024, fasta);
    return (Box::new(reader), pb)
}

fn open_file(filename: String) -> Box<dyn Read>
{
    let file = match File::open(&filename) {
        Err(why) => panic!("Couldn't open {}: {}", filename, why.to_string()),
        Ok(file) => file,
    };

    let file = BufReader::with_capacity(64 * 1024 * 1024, file);

    let fasta: Box<dyn Read> = if filename.ends_with("gz") {
        Box::new(flate2::read::GzDecoder::new(file))
    } else if filename.ends_with("snappy") || filename.ends_with("sz") {
        Box::new(snap::read::FrameDecoder::new(file))
    } else {
        Box::new(file)
    };

    let reader = BufReader::with_capacity(32 * 1024 * 1024, fasta);
    return Box::new(reader)
}

impl Iterator for Sequences {
    type Item = Sequence;

    fn next(&mut self) -> Option<Sequence> {
        let seq: Sequence = match bincode::deserialize_from(&mut self.reader) {
            Ok(x)   => x,
            Err(_)  => return None
        };
        Some(seq)
    }
}

impl Sequences {
    pub fn new(filename: String) -> Sequences {
        return Sequences { reader: open_file(filename) }
    }
}