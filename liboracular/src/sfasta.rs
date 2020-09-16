use std::io::prelude::*;

use std::convert::From;

use std::fs::{metadata, File};
use std::io::{BufRead, BufReader, BufWriter, Read, SeekFrom, Write};
use std::path::Path;
use std::time::Instant;

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use crate::io;

// TODO: Set a const for BufReader buffer size
//       Make it a global const, but also maybe make it configurable?
//       Reason being that network FS will benefit from larger buffers
// TODO: Also make BufWriter bufsize global, but ok to leave larger.

/// Represents an entry from an SFASTA file
#[derive(PartialEq, Serialize, Deserialize, Debug)]
pub struct Entry {
    pub id: String,
    pub seq: Vec<u8>,
}

/// SFASTA files stored on disk are bincoded, with the sequence being
/// compressed. Decompression does not occur unless EntryCompressed is converted
/// to an Entry. This allows faster searching of SFASTA files without spending
/// CPU cycles on decompression prematurely.
#[derive(PartialEq, Serialize, Deserialize)]
pub struct EntryCompressed {
    pub id: String,
    pub compressed_seq: Vec<u8>,
}

/// Destroys EntryCompressed struct and returns owned id as String
impl EntryCompressed {
    pub fn take_id(self) -> String {
        self.id
    }
}

/// Iterator to return io::Sequences
pub struct Sequences {
    reader: Box<dyn Read + Send>,
}

// TODO: Should be a TryFrom really...
/// Convert an EntryCompressed to Entry.
/// Automatically performs the decompression.
impl From<EntryCompressed> for Entry {
    fn from(item: EntryCompressed) -> Self {
        let len = item.compressed_seq.len();
        let mut seq_reader = snap::read::FrameDecoder::new(&item.compressed_seq[..]);
        let mut seq: Vec<u8> = Vec::with_capacity(len * 2);
        seq_reader
            .read_to_end(&mut seq)
            .expect("Unable to read compressed sequence");
        Entry { id: item.id, seq }
    }
}

// TODO: Should be a TryFrom really...
/// Converts an Entry into EntryCompressed. Performs the compression
/// automatically.
impl From<Entry> for EntryCompressed {
    fn from(item: Entry) -> Self {
        let compressed_seq = Vec::with_capacity(item.seq.len());
        let mut writer = snap::write::FrameEncoder::new(compressed_seq);
        writer
            .write_all(&item.seq)
            .expect("Unable to write to vector...");
        writer.flush().expect("Unable to flush");
        let compressed_seq = writer.into_inner().unwrap();
        EntryCompressed {
            id: item.id,
            compressed_seq,
        }
    }
}

/// Converts an SFASTA::Entry into io::Sequence for further processing
/// Really an identity function...
impl From<Entry> for io::Sequence {
    fn from(item: Entry) -> Self {
        let len = item.seq.len();
        io::Sequence {
            id: item.id,
            seq: item.seq,
            location: 0,
            end: len,
        }
    }
}

impl Sequences {
    /// Given a filename, returns a Sequences variable.
    /// Can be used as an iterator.
    pub fn new(filename: String) -> Sequences {
        Sequences {
            reader: open_file(filename),
        }
    }
}

// TODO: Create the option to pass back lowercase stuff too..
// Maybe for repeat masking and such? Right now it's all uppercase.
impl Iterator for Sequences {
    type Item = io::Sequence;

    /// Get the next SFASTA entry as io::Sequence type
    fn next(&mut self) -> Option<io::Sequence> {
        let ec: EntryCompressed = match bincode::deserialize_from(&mut self.reader) {
            Ok(x) => x,
            Err(_) => return None, // panic!("Error at SFASTA::Sequences::next: {}", y)
        };

        // Have to convert from EntryCompressed to Entry, this handles that middle
        // conversion.
        let middle: Entry = ec.into();
        let mut seq: io::Sequence = middle.into();
        seq.make_uppercase();

        Some(seq)
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
    writer
        .write_all(&sequence)
        .expect("Unable to write to vector...");
    writer.flush().expect("Unable to flush");
    compressed = writer.into_inner().unwrap();

    EntryCompressed {
        id,
        compressed_seq: compressed,
    }
}

/// Converts a FASTA file to an SFASTA file...
pub fn convert_fasta_file(filename: String, output: String)
// TODO: Add progress bar option
//
// Convert file to bincode/snappy for faster processing
// Stores accession/taxon information inside the Sequence struct
{
    let mut total_counts: usize = 0;

    let output_filename = check_extension(output);

    let out_file = File::create(output_filename.clone()).expect("Unable to write to file");
    let mut out_fh = BufWriter::with_capacity(4 * 1024 * 1024, out_file);

    let mut buffer: Vec<u8> = Vec::with_capacity(1024);
    let mut id: String = String::from("INVALID_ID_FIRST_ENTRY_YOU_SHOULD_NOT_SEE_THIS");
    let mut seqbuffer: Vec<u8> = Vec::with_capacity(32 * 1024 * 1024);
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

    let mut reader = BufReader::with_capacity(512 * 1024, fasta);
    while let Ok(bytes_read) = reader.read_until(b'\n', &mut buffer) {
        if bytes_read == 0 {
            // No more reads, thus no more data...
            // Write out the final piece...
            let entry = generate_sfasta_entry(id, seqbuffer[..seqlen].to_vec());
            bincode::serialize_into(&mut out_fh, &entry)
                .expect("Unable to write to bincode output");

            break;
        }

        match buffer[0] {
            // 62 is a > meaning we have a new sequence id.
            62 => {
                total_counts += 1;
                // Write out entry...
                if seqlen > 0 {
                    // Ignore first one..
                    let entry = generate_sfasta_entry(id, seqbuffer[..seqlen].to_vec());
                    bincode::serialize_into(&mut out_fh, &entry)
                        .expect("Unable to write to bincode output");

                    seqbuffer.clear();
                    seqlen = 0;
                }

                let slice_end = bytes_read.saturating_sub(1);
                id = String::from_utf8(buffer[1..slice_end].to_vec())
                    .expect("Invalid UTF-8 encoding...");
                id = id.split(' ').next().unwrap().trim().to_string();
            }
            _ => {
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

    let mut reader = BufReader::with_capacity(512 * 1024, file);

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
pub fn get_headers_from_sfasta(filename: String) -> Vec<String> {
    let file = match File::open(&filename) {
        Err(why) => panic!("Couldn't open {}: {}", filename, why.to_string()),
        Ok(file) => file,
    };

    let mut reader = BufReader::with_capacity(512 * 1024, file);

    let mut ids: Vec<String> = Vec::with_capacity(2048);

    while let Ok(entry) = bincode::deserialize_from::<_, EntryCompressed>(&mut reader) {
        ids.push(entry.id);
    }

    ids
}

/// Get all IDs from an SFASTA file
/// Really a debugging function...
pub fn test_sfasta(filename: String) {
    let file = match File::open(&filename) {
        Err(why) => panic!("Couldn't open {}: {}", filename, why.to_string()),
        Ok(file) => file,
    };

    let mut reader = BufReader::with_capacity(512 * 1024, file);

    let mut seqnum: usize = 0;

    loop {
        seqnum += 1;
        match bincode::deserialize_from::<_, EntryCompressed>(&mut reader) {
            Ok(_) => println!("OK SEQ: {}", seqnum),
            Err(x) => panic!("Found error: {}", x),
        };
    }
}

/// Checks that the file extension ends in .sfasta or adds it if necessary
#[inline(always)]
fn check_extension(filename: String) -> String {
    if !filename.ends_with(".sfasta") {
        format!("{}.sfasta", filename)
    } else {
        filename
    }
}

/// Opens an SFASTA file and returns a Box<dyn Read> type
fn open_file(filename: String) -> Box<dyn Read + Send> {
    let filename = check_extension(filename);

    let file = match File::open(&filename) {
        Err(_) => panic!("Couldn't open {}", filename),
        Ok(file) => file,
    };

    let reader = BufReader::with_capacity(512 * 1024, file);

    Box::new(reader)
}

/// Indexes an SFASTA file
pub fn index(filename: &str) -> String {
    // TODO: Run a sanity check on the file first... Make sure it's valid
    // sfasta

    let filesize = metadata(&filename).expect("Unable to open file").len();
    let starting_size = std::cmp::max((filesize / 1000) as usize, 1024);
    println!("Starting Size: {}", starting_size);

    //    let mut idx: HashMap<String, u64, RandomXxHashBuilder64> =
    // Default::default();    idx.reserve(starting_size);

    let mut ids = Vec::with_capacity(starting_size);
    let mut locations = Vec::with_capacity(starting_size);

    let fh = match File::open(&filename) {
        Err(why) => panic!("Couldn't open {}: {}", filename, why.to_string()),
        Ok(file) => file,
    };

    let mut fh = BufReader::with_capacity(512 * 1024, fh);
    let mut pos = fh
        .seek(SeekFrom::Current(0))
        .expect("Unable to work with seek API");
    let mut i = 0;
    let mut now = Instant::now();
    //    let mut bump = Bump::new();
    //    let mut maxalloc: usize = 0;

    while let Ok(entry) = bincode::deserialize_from::<_, EntryCompressed>(&mut fh) {
        i += 1;
        if i % 100_000 == 0 {
            println!("100k at {} ms.", now.elapsed().as_millis()); //Maxalloc {} bytes", now.elapsed().as_secs(), maxalloc);
            println!(
                "{}/{} {}",
                pos,
                filesize,
                (pos as f32 / filesize as f32) as f32
            );
            now = Instant::now();
        }

        ids.push(entry.take_id());
        locations.push(pos);
        //        idx.insert(entry.id.clone(), pos);
        pos = fh
            .seek(SeekFrom::Current(0))
            .expect("Unable to work with seek API");
        //        maxalloc = std::cmp::max(maxalloc, bump.allocated_bytes());
        //        bump.reset();
    }
    println!("Finished with {} steps", i);

    let idx: HashMap<String, u64> = ids.into_iter().zip(locations).collect();

    let filenamepath = Path::new(&filename);
    let filename = Path::new(filenamepath.file_name().unwrap())
        .file_stem()
        .unwrap()
        .to_str()
        .unwrap()
        .to_owned()
        + ".sfai";
    let output_filename =
        filenamepath.parent().unwrap().to_str().unwrap().to_owned() + "/" + &filename;

    println!("Saving to file: {}", output_filename);

    let out_file = snap::write::FrameEncoder::new(
        File::create(output_filename.clone()).expect("Unable to write to file"),
    );
    let mut out_fh = BufWriter::with_capacity(4 * 1024 * 1024, out_file);
    bincode::serialize_into(&mut out_fh, &idx).expect("Unable to write index");

    output_filename
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::prelude::*;

    #[test]
    pub fn convert_fasta_to_sfasta_and_index() {
        let input_filename = "test_data/test_multiple.fna";
        let output_filename = "test_data/test_sfasta_convert_and_index.sfasta";

        convert_fasta_file(input_filename.to_string(), output_filename.to_string());

        let idx_filename = index(output_filename);
        assert!(idx_filename == "test_data/test_sfasta_convert_and_index.sfai");
    }
}
