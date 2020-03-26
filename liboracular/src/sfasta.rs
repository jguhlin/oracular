use std::io::prelude::*;
use std::io::BufWriter;

use std::sync::{Arc, RwLock};

use std::fs::File;
use std::io::{BufReader, Read, BufRead};

use serde::{Serialize, Deserialize};

#[derive(PartialEq, Serialize, Deserialize)]
pub struct Entry {
    pub id:  String,
    pub seq: Vec<u8>,
}

#[derive(PartialEq, Serialize, Deserialize)]
pub struct EntryCompressed {
    pub id:  String,
    pub compressed_seq: Vec<u8>,
}

// sfasta is:
// bincode encoded
//   fasta ID
//   snappy compressed sequence
// TODO: Add some form of indexing?

// TODO: Make into trait and impl and From and To stuff...
fn generate_sfasta_entry(id: String, sequence: Vec<u8>) -> EntryCompressed {
    let mut compressed: Vec<u8> = Vec::with_capacity(sequence.len());
    let mut writer = snap::write::FrameEncoder::new(compressed);
    writer.write_all(&sequence).expect("Unable to write to vector...");
    writer.flush().expect("Unable to flush");
    compressed = writer.into_inner().unwrap();

    EntryCompressed { id, compressed_seq: compressed }
}

pub fn convert_fasta_file(filename: String, output: String,)
// TODO: Add progress bar option
// Convert file to bincode/snappy for faster processing
// Stores accession/taxon information inside the Sequence struct
{

    let output_filename;
    if !output.ends_with(".sfasta") {
        output_filename = format!("{}.sfasta", output);
    } else {
        output_filename = output;
    }

    let out_file = File::create(output_filename).expect("Unable to write to file");
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

    let mut reader = BufReader::with_capacity(64 * 1024 * 1024, fasta);
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
}

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