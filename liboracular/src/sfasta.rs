use std::io::prelude::*;
use std::io::BufWriter;

use std::sync::{Arc, RwLock};

use std::convert::From;

use std::fs::File;
use std::io::{BufReader, Read, BufRead};

use serde::{Serialize, Deserialize};

use crate::io;

/// Represents an entry from an SFASTA file
#[derive(PartialEq, Serialize, Deserialize)]
pub struct Entry {
    pub id:  String,
    pub seq: Vec<u8>,
}

/// SFASTA files stored on disk are bincoded, with the sequence being compressed.
/// Decompression does not occur unless EntryCompressed is converted to an Entry.
/// This allows faster searching of SFASTA files without spending CPU cycles on
/// decompression prematurely.
#[derive(PartialEq, Serialize, Deserialize)]
pub struct EntryCompressed {
    pub id:  String,
    pub compressed_seq: Vec<u8>,
}

/// Iterator to return io::Sequences
pub struct Sequences {
    reader: Box<dyn Read>,
}

// TODO: Should be a TryFrom really...
/// Convert an EntryCompressed to Entry.
/// Automatically performs the decompression.
impl From<EntryCompressed> for Entry {
    fn from(item: EntryCompressed) -> Self {
        let len = item.compressed_seq.len();
        let mut seq_reader = snap::read::FrameDecoder::new(&item.compressed_seq[..]);
        let mut seq: Vec<u8> = Vec::with_capacity(len * 2);
        seq_reader.read_to_end(&mut seq).expect("Unable to read compressed sequence");
        Entry { id: item.id, seq }
    }
}

// TODO: Should be a TryFrom really...
/// Converts an Entry into EntryCompressed. Performs the compression automatically.
impl From<Entry> for EntryCompressed {
    fn from(item: Entry) -> Self {
        let compressed_seq = Vec::with_capacity(item.seq.len());
        let mut writer = snap::write::FrameEncoder::new(compressed_seq);
        writer.write_all(&item.seq).expect("Unable to write to vector...");
        writer.flush().expect("Unable to flush");
        let compressed_seq = writer.into_inner().unwrap();
        EntryCompressed{ id: item.id, compressed_seq }
    }
}

/// Converts an SFASTA::Entry into io::Sequence for further processing
/// Really an identity function...
impl From<Entry> for io::Sequence {
    fn from(item: Entry) -> Self {
        io::Sequence { id: item.id, seq: item.seq }
    }
}

impl Sequences {
    /// Given a filename, returns a Sequences variable.
    /// Can be used as an iterator.
    pub fn new(filename: String) -> Sequences {
        return Sequences { reader: open_file(filename) }
    }
}

impl Iterator for Sequences {
    type Item = io::Sequence;

    /// Get the next SFASTA entry as io::Sequence type
    fn next(&mut self) -> Option<io::Sequence> {
        let ec: EntryCompressed = match bincode::deserialize_from(&mut self.reader) {
            Ok(x)   => x,
            Err(y)  => return None // panic!("Error at SFASTA::Sequences::next: {}", y)
        };

        // Have to convert from EntryCompressed to Entry, this handles that middle
        // conversion.
        let middle: Entry = ec.into();

        Some(middle.into())
    }
}

// sfasta is:
// bincode encoded
//   fasta ID
//   snappy compressed sequence
// TODO: Add some form of indexing?

// TODO: Remove this code since we have the trait now...
/// Should remove...
fn generate_sfasta_entry(id: String, sequence: Vec<u8>) -> EntryCompressed {
    let mut compressed: Vec<u8> = Vec::with_capacity(sequence.len());
    let mut writer = snap::write::FrameEncoder::new(compressed);
    writer.write_all(&sequence).expect("Unable to write to vector...");
    writer.flush().expect("Unable to flush");
    compressed = writer.into_inner().unwrap();

    EntryCompressed { id, compressed_seq: compressed }
}

/// Converts a FASTA file to an SFASTA file...
pub fn convert_fasta_file(filename: String, output: String,)
// TODO: Add progress bar option
//
// Convert file to bincode/snappy for faster processing
// Stores accession/taxon information inside the Sequence struct
{
    let mut total_counts: usize = 0;

    let output_filename = check_extension(output);

    let out_file = File::create(output_filename.clone()).expect("Unable to write to file");
    let mut out_fh = BufWriter::with_capacity(64 * 1024 * 1024, out_file);

    let mut buffer: Vec<u8> = Vec::with_capacity(1024);
    let mut id: String = String::from("INVALID_ID_FIRST_ENTRY_YOU_SHOULD_NOT_SEE_THIS");
    let mut seqbuffer: Vec<u8> = Vec::with_capacity(64 * 1024 * 1024);
    let mut seqlen: usize = 0;

    let file = match File::open(&filename) {
        Err(why) => panic!("Couldn't open {}: {}", filename, why.to_string()),
        Ok(file) => file,
    };

    let fasta: Box<dyn Read> = if filename.ends_with("gz") {
        Box::new(flate2::read::GzDecoder::new(file))
    } else if filename.ends_with("snappy") || filename.ends_with("sz") {
        Box::new(snap::read::FrameDecoder::new(file))
    } else {
        Box::new(file)
    };

    let mut reader = BufReader::with_capacity(32 * 1024 * 1024, fasta);
    while let Ok(bytes_read) = reader.read_until(b'\n', &mut buffer) {

        if bytes_read == 0 { 
            // No more reads, thus no more data...
            // Write out the final piece...
            let entry = generate_sfasta_entry(id, seqbuffer[..seqlen].to_vec());
            bincode::serialize_into(&mut out_fh, &entry).expect("Unable to write to bincode output");

            break;
        }

        match buffer[0] {
            // 62 is a > meaning we have a new sequence id.
            62 => {
                total_counts += 1;
                // Write out entry...
                if seqlen > 0 { // Ignore first one..
                    let entry = generate_sfasta_entry(id, seqbuffer[..seqlen].to_vec());
                    bincode::serialize_into(&mut out_fh, &entry).expect("Unable to write to bincode output");

                    seqbuffer.clear();
                    seqlen = 0;
                }

                let slice_end = bytes_read.saturating_sub(1);
                id = String::from_utf8(buffer[1..slice_end].to_vec()).expect("Invalid UTF-8 encoding...");
            },
            _  => {
                let slice_end = bytes_read.saturating_sub(1);
                seqbuffer.extend_from_slice(&buffer[0..slice_end]);
                seqlen = seqlen.saturating_add(slice_end);
            }
        }

        buffer.clear();
    }

    drop(out_fh);

    let file = match File::open(&output_filename) {
        Err(why) => panic!("Couldn't open {}: {}", filename, why.to_string()),
        Ok(file) => file,
    };

    let mut reader = BufReader::with_capacity(32 * 1024 * 1024, file);

    let mut ids: Vec<String> = Vec::with_capacity(2048);

    let mut saved_count: usize = 0;

    let mut last_entry: String = "".to_string();

    while let Ok(entry) = bincode::deserialize_from::<_, EntryCompressed>(&mut reader) {
        saved_count += 1;
        last_entry = entry.id;
    }

    if saved_count != total_counts {
        println!("Error: Not the right amount!");
        println!("Last Entry: {}", last_entry);
        panic!("Error converting file");
    }
}

/// Get all IDs from an SFASTA file
/// Really a debugging function...
pub fn get_headers_from_sfasta(filename: String) -> Vec<String>
{
    let file = match File::open(&filename) {
        Err(why) => panic!("Couldn't open {}: {}", filename, why.to_string()),
        Ok(file) => file,
    };

    let mut reader = BufReader::with_capacity(32 * 1024 * 1024, file);

    let mut ids: Vec<String> = Vec::with_capacity(2048);

    while let Ok(entry) = bincode::deserialize_from::<_, EntryCompressed>(&mut reader) {
        ids.push(entry.id);
    }

    return ids
}

/// Get all IDs from an SFASTA file
/// Really a debugging function...
pub fn test_sfasta(filename: String)
{
    let file = match File::open(&filename) {
        Err(why) => panic!("Couldn't open {}: {}", filename, why.to_string()),
        Ok(file) => file,
    };

    let mut reader = BufReader::with_capacity(32 * 1024 * 1024, file);

    loop {
        match bincode::deserialize_from::<_, EntryCompressed>(&mut reader) {
            Ok(entry) => (),
            Err(x)    => panic!("Found error: {}", x),
        };
    }
}


/// Checks that the file extension ends in .sfasta or adds it if necessary
#[inline(always)]
fn check_extension(filename: String) -> String {
    let outfname;
    if !filename.ends_with(".sfasta") {
        outfname = format!("{}.sfasta", filename);
    } else {
        outfname = filename;
    }
    outfname
}

/// Opens an SFASTA file and returns a Box<dyn Read> type
fn open_file(filename: String) -> Box<dyn Read> {

    let filename = check_extension(filename);

    let file = match File::open(&filename) {
        Err(why) => panic!("Couldn't open {}: {}", filename, why.to_string()),
        Ok(file) => file,
    };

    let file = BufReader::with_capacity(32 * 1024 * 1024, file);
    let reader = BufReader::with_capacity(32 * 1024 * 1024, file);

    return Box::new(reader)
}

/*
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

    let mut reader = BufReader::with_capacity(32 * 1024 * 1024, fasta);
    return (Box::new(reader), pb)
}

fn snappy_output(filename: String) -> Box<dyn Write> {
    let buffer = BufWriter::with_capacity(64 * 1024 * 1024, 
        File::create(filename).expect("Unable to write to file!"));
    Box::new(BufWriter::with_capacity(16 * 1024 * 1024, snap::write::FrameEncoder::new(buffer)))
}*/