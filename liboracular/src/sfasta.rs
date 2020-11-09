use hashbrown::HashMap;
use std::convert::From;
use std::convert::TryFrom;
use std::fs::{metadata, File};
use std::io::prelude::*;
use std::io::{BufRead, BufReader, BufWriter, Read, SeekFrom, Write};
use std::path::Path;
use std::time::Instant;

use crate::utils::generic_open_file;
use crate::fasta;

use rand::prelude::*;
use rand_chacha::ChaCha20Rng;

use serde::{Deserialize, Serialize};

use crate::io;

// SuperTrait
pub trait ReadAndSeek: Read + Seek + Send {}
impl<T: Read + Seek + Send> ReadAndSeek for T {}

// TODO: Spin this out as a separate library...
// TODO: Set a const for BufReader buffer size
//       Make it a global const, but also maybe make it configurable?
//       Reason being that network FS will benefit from larger buffers
// TODO: Also make BufWriter bufsize global, but ok to leave larger.
#[derive(PartialEq, Serialize, Deserialize, Debug, Clone)]
pub enum CompressionType {
    ZSTD,
    SNAPPY,
    GZIP, // Please don't use this -- IMPLEMENT
    NAF,  // Not yet supported -- IMPLEMENT
    NONE, // No Compression -- IMPLEMENT
}

#[derive(PartialEq)]
pub enum SeqMode {
    Linear,
    Random,
}

#[derive(PartialEq, Serialize, Deserialize, Debug, Clone)]
pub struct Header {
    pub id: Option<String>,
    pub comment: Option<String>,
    pub citation: Option<String>,
    #[serde(with = "serde_bytes")]
    pub dict: Option<Vec<u8>>,
}

/// Represents an entry from an SFASTA file (extension of .sfasta)
#[derive(PartialEq, Serialize, Deserialize, Debug)]
pub struct Entry {
    pub id: String,
    #[serde(with = "serde_bytes")]
    pub seq: Vec<u8>,
    pub comment: Option<String>,
    pub len: u64,
}

impl Entry {
    pub fn compress(
        self,
        compression_type: CompressionType,
        compression_level: i32,
        dict: &Option<Vec<u8>>,
    ) -> EntryCompressed {
        let Entry {
            id,
            seq,
            comment,
            len,
        } = self;

        let compressed_seq = Vec::with_capacity(self.len as usize);
        // SNAPPY Line
        //let mut writer = snap::write::FrameEncoder::new(compressed_seq);

        // ZSTD Line
        let mut writer = if dict.is_some() {
                zstd::stream::Encoder::with_dictionary(
                    compressed_seq,
                    compression_level,
                    &dict.as_ref().unwrap(),
                )
                .unwrap()
            } else {
                zstd::stream::Encoder::new(compressed_seq, compression_level).unwrap()
            };
        writer
            .write_all(&seq)
            .expect("Unable to write to vector...");
        writer.flush().expect("Unable to flush");
        let compressed_seq = writer.finish().unwrap();

        EntryCompressed {
            id,
            compression_type,
            compressed_seq,
            comment,
            len,
        }
    }
}

/// SFASTA files stored on disk are bincoded, with the sequence being
/// compressed. Decompression does not occur unless EntryCompressed is converted
/// to an Entry. This allows faster searching of SFASTA files without spending
/// CPU cycles on decompression prematurely.
#[derive(PartialEq, Serialize, Deserialize, Clone)]
pub struct EntryCompressed {
    pub id: String,
    pub compression_type: CompressionType,

    #[serde(with = "serde_bytes")]
    pub compressed_seq: Vec<u8>,
    pub comment: Option<String>,
    pub len: u64,
}

/// Destroys EntryCompressed struct and returns owned id as String
impl EntryCompressed {
    pub fn take_id(self) -> String {
        self.id
    }

    pub fn decompress(self, dict: &Option<Vec<u8>>) -> Entry {
        let EntryCompressed {
            id,
            compression_type,
            compressed_seq,
            comment,
            len,
        } = self;
        // SNAPPY Compatability line
        //let mut seq_reader = snap::read::FrameDecoder::new(&item.compressed_seq[..]);

        let mut seq: Vec<u8> = Vec::with_capacity(len as usize);

        // ZSTD Lines
        if dict.is_some() {
            let mut seq_reader = zstd::stream::read::Decoder::with_dictionary(
                &compressed_seq[..],
                dict.as_ref().unwrap(),
            )
            .unwrap();
            seq_reader
                .read_to_end(&mut seq)
                .expect("Unable to read compressed sequence");
        } else {
            let mut seq_reader = zstd::stream::read::Decoder::new(&compressed_seq[..]).unwrap();
            seq_reader
                .read_to_end(&mut seq)
                .expect("Unable to read compressed sequence");
        }

        Entry {
            id,
            seq,
            comment,
            len,
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

/// Iterator to return sfasta::EntryCompressed
pub struct CompressedSequences {
    pub header: Header,
    reader: Box<dyn ReadAndSeek + Send>,
    pub idx: Option<(Vec<String>, Vec<u64>)>,
    access: SeqMode,
    random_list: Option<Vec<u64>>,
}

/// Iterator to return io::Sequences
pub struct Sequences {
    pub header: Header,
    reader: Box<dyn ReadAndSeek + Send>,
    pub idx: Option<(Vec<String>, Vec<u64>)>,
    access: SeqMode,
    random_list: Option<Vec<u64>>,
}

// TODO: We need to cache DecoderDictionary at some point...
impl Sequences {
    /// Given a filename, returns a Sequences variable.
    /// Can be used as an iterator.
    pub fn new(filename: &str) -> Sequences {
        let (mut reader, idx) = open_file(filename);
        let header: Header = match bincode::deserialize_from(&mut reader) {
            Ok(x) => x,
            Err(_) => panic!("Header missing or malformed in SFASTA file"),
        };

        Sequences {
            header,
            reader,
            idx,
            access: SeqMode::Linear,
            random_list: None,
        }
    }

    pub fn set_mode(&mut self, access: SeqMode) {
        self.access = access;

        if self.access == SeqMode::Random {
            assert!(self.idx.is_some());
            let idx = self.idx.as_ref().unwrap();
            let mut locs: Vec<u64> = idx.1.clone();
            let mut rng = ChaCha20Rng::seed_from_u64(42);
            locs.shuffle(&mut rng);
            self.random_list = Some(locs);
        } else {
            self.random_list = None;
        }
    }

    // Convert to iterator that only returns EntryCompressed...
    pub fn into_compressed_sequences(self) -> CompressedSequences {
        let Sequences {
            header,
            reader,
            idx,
            access,
            random_list,
        } = self;

        CompressedSequences {
            header,
            reader,
            idx,
            access,
            random_list,
        }
    }
}

// TODO: Create the option to pass back lowercase stuff too..
// Maybe for repeat masking and such? Right now it's all uppercase.
impl Iterator for Sequences {
    type Item = io::Sequence;

    /// Get the next SFASTA entry as io::Sequence type
    fn next(&mut self) -> Option<io::Sequence> {
        if self.access == SeqMode::Random {
            assert!(self.random_list.is_some());
            let rl = self.random_list.as_mut().unwrap();
            let next_loc = rl.pop();
            next_loc?;
            self.reader
                .seek(SeekFrom::Start(next_loc.unwrap()))
                .expect("Unable to work with seek API");
        }

        let ec: EntryCompressed = match bincode::deserialize_from(&mut self.reader) {
            Ok(x) => x,
            Err(_) => return None, // panic!("Error at SFASTA::Sequences::next: {}", y)
        };

        // Have to convert from EntryCompressed to Entry, this handles that middle
        // conversion.
        let middle: Entry = ec.decompress(&self.header.dict);
        let mut seq: io::Sequence = middle.into();
        seq.make_uppercase();

        Some(seq)
    }
}

impl Iterator for CompressedSequences {
    type Item = EntryCompressed;

    /// Get the next SFASTA entry as an EntryCompressed struct
    fn next(&mut self) -> Option<EntryCompressed> {
        if self.access == SeqMode::Random {
            assert!(self.random_list.is_some());
            let rl = self.random_list.as_mut().unwrap();
            let next_loc = rl.pop();
            next_loc?;
            self.reader
                .seek(SeekFrom::Start(next_loc.unwrap()))
                .expect("Unable to work with seek API");
        }

        let ec: EntryCompressed = match bincode::deserialize_from(&mut self.reader) {
            Ok(x) => x,
            Err(_) => return None, // panic!("Error at SFASTA::Sequences::next: {}", y)
        };

        // Have to convert from EntryCompressed to Entry, this handles that middle
        // conversion.
        Some(ec)
    }
}

// sfasta is:
// bincode encoded
//   fasta ID
//   zstd compressed sequence

// TODO: Remove this code since we have the trait now...
/// Should remove...
fn generate_sfasta_compressed_entry(
    id: String,
    comment: Option<String>,
    seq: Vec<u8>,
    compression_type: CompressionType,
    compression_level: i32,
    dict: &Option<Vec<u8>>,
) -> EntryCompressed {
    let len: u64 =
        u64::try_from(seq.len()).expect("Unlikely as it is, sequence length exceeds u64::MAX...");
    let entry = Entry {
        id,
        seq,
        comment,
        len,
    };

    entry.compress(compression_type, compression_level, dict)
}

/// Opens an SFASTA file and an index and returns a Box<dyn Read>,
/// HashMap<String, usize> type
pub fn open_file(filename: &str) -> (Box<dyn ReadAndSeek + Send>, Option<(Vec<String>, Vec<u64>)>) {
    let filename = check_extension(filename);

    let file = match File::open(Path::new(&filename)) {
        Err(_) => panic!("Couldn't open {}", filename),
        Ok(file) => file,
    };

    let reader = BufReader::with_capacity(512 * 1024, file);

    (Box::new(reader), load_index(&filename))
}

/// Build a ZSTD dictionary
pub fn build_zstd_dict(filename: &str) -> Option<Vec<u8>> {
    let (_, _, fasta) = generic_open_file(filename);
    let mut reader = BufReader::with_capacity(512 * 1024, fasta);

    let mut buffer: Vec<u8> = Vec::with_capacity(256);

    let mut sample_data: Vec<u8> = Vec::with_capacity(64 * 1024 * 1024);
    let mut sample_sizes: Vec<usize> = Vec::with_capacity(1024);
    let mut maxsize: usize = 0;
    let mut total_len: usize = 0;
    let mut cur_len: usize = 0;
    let mut first: bool = true;

    while let Ok(bytes_read) = reader.read_until(b'\n', &mut buffer) {
        if bytes_read == 0 || total_len >= 64 * 1024 * 1024 {
            maxsize = std::cmp::max(maxsize, cur_len);
            sample_sizes.push(cur_len);
            break;
        }

        match buffer[0] {
            // 62 is a > meaning we have a new sequence id OR this is the first entry...
            62 => {
                if first {
                    first = false;
                } else {
                    maxsize = std::cmp::max(maxsize, cur_len);
                    sample_sizes.push(cur_len);
                    cur_len = 0;
                }
            }
            _ => {
                let slice_end = bytes_read.saturating_sub(1);
                sample_data.extend_from_slice(&buffer[0..slice_end]);
                cur_len += slice_end;
                total_len += slice_end;
            }
        }

        buffer.clear();
    }

    assert!(!sample_sizes.is_empty());

    if total_len <= 1024 * 256 || sample_sizes.len() <= 7 {
        None
    } else {
        match zstd::dict::from_continuous(&sample_data, &sample_sizes, maxsize) {
            Ok(x) => Some(x),
            Err(_) => None,
        }
    }
}

use std::thread;
use std::thread::{park, JoinHandle};
use std::sync::atomic::Ordering;
use std::sync::atomic::{AtomicBool, AtomicUsize};
use std::sync::Arc;
use crossbeam::queue::ArrayQueue;
use crossbeam::utils::Backoff;

/// Converts a FASTA file to an SFASTA file...
pub fn convert_fasta_file(filename: &str, output: &str)
// TODO: Make multithreaded for very large datasets (>= 1Gbp or 5Gbp or something)
// TODO: Add progress bar option ... or not..
//
// Convert file to bincode/zstd for faster processing
// Stores accession/taxon information inside the Sequence struct
{
    // let dict = build_zstd_dict(filename);
    let dict = None;

    let output_filename = check_extension(output);

    let out_file = File::create(output_filename.clone()).expect("Unable to write to file");
    let mut out_fh = BufWriter::with_capacity(1024 * 1024, out_file);

    let (filesize, _, _) = generic_open_file(filename);

    let header = Header {
        dict,
        citation: None,
        comment: None,
        id: Some(filename.to_string()),
    };

    bincode::serialize_into(&mut out_fh, &header).expect("Unable to write to bincode output");

    let starting_size = std::cmp::max((filesize / 500) as usize, 1024);

    let mut ids = Vec::with_capacity(starting_size);
    let mut locations = Vec::with_capacity(starting_size);

    let mut pos = out_fh
        .seek(SeekFrom::Current(0))
        .expect("Unable to work with seek API");

    let fasta = fasta::Fasta::new(filename);

    let thread_count = 64;
    let queue_size = 1024;

    // multi-threading...
    let shutdown = Arc::new(AtomicBool::new(false));
    let total_entries = Arc::new(AtomicUsize::new(0));
    let compressed_entries = Arc::new(AtomicUsize::new(0));
    let written_entries = Arc::new(AtomicUsize::new(0));

    let queue: Arc<ArrayQueue<fasta::Sequence>> = Arc::new(ArrayQueue::new(queue_size));
    let output_queue: Arc<ArrayQueue<EntryCompressed>> = Arc::new(ArrayQueue::new(queue_size));

    let mut worker_handles = Vec::new();

    for _ in 0..thread_count {
        let q = Arc::clone(&queue);
        let oq = Arc::clone(&output_queue);
        let header_copy = header.clone();
        let shutdown_copy = Arc::clone(&shutdown);
        let te = Arc::clone(&total_entries);
        let ce = Arc::clone(&compressed_entries);

        let handle = thread::spawn(move || {
            let shutdown = shutdown_copy;
            let backoff = Backoff::new();
            let header = header_copy;
            let mut result;
            loop {
                result = q.pop();

                match result {
                    None => {
                        backoff.snooze();
                        if shutdown.load(Ordering::Relaxed) && ce.load(Ordering::Relaxed) == te.load(Ordering::Relaxed) {
                            return;
                        }
                    },
                    Some(x) => {
                        // let x = result.unwrap();
                        let mut entry: EntryCompressed = generate_sfasta_compressed_entry(
                            x.id.clone(),
                            None,
                            x.seq.to_vec(),
                            CompressionType::ZSTD,
                            3,
                            &header.dict,
                        );

                        while let Err(x) = oq.push(entry) {
                            entry = x;
                            park(); // Queue is full, park the thread...
                        }
                        ce.fetch_add(1, Ordering::SeqCst);
                    }
                }
            }
        });

        worker_handles.push(handle);
    }

    let oq = Arc::clone(&output_queue);
    let shutdown_copy = Arc::clone(&shutdown);
    let q = Arc::clone(&queue);
    let te = Arc::clone(&total_entries);

    let output_thread = thread::spawn(move || {
        let shutdown = shutdown_copy;
        let output_queue = oq;
        let backoff = Backoff::new();

        let mut result;
        loop {
            result = output_queue.pop();
            match result {
                None => {
                    // Unpark all other threads..
                    for i in &worker_handles {
                        i.thread().unpark();
                    }
                    backoff.snooze();
                    if (written_entries.load(Ordering::Relaxed) == te.load(Ordering::Relaxed))
                        && shutdown.load(Ordering::Relaxed) {
                        drop(out_fh);
                        create_index(&output_filename, ids, locations);
                        return;
                    }
                },
                Some(cs) => {
                    ids.push(cs.id.clone());
                    locations.push(pos);
                    bincode::serialize_into(&mut out_fh, &cs)
                        .expect("Unable to write to bincode output");
                    pos = out_fh
                        .seek(SeekFrom::Current(0))
                        .expect("Unable to work with seek API");
                    written_entries.fetch_add(1, Ordering::SeqCst);
                }
            }
        }
    });

    let backoff = Backoff::new();
    fasta.for_each(|x| {
        let mut item;
        item = x;
        while let Err(x) = queue.push(item) {
            item = x;
            backoff.snooze();
        }
        total_entries.fetch_add(1, Ordering::SeqCst);

    });

    while queue.len() > 0 || output_queue.len() > 0 {
        backoff.snooze();
    }

    shutdown.store(true, Ordering::SeqCst);

    output_thread.join().expect("Unable to join the output thread back...");
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
fn check_extension(filename: &str) -> String {
    if !filename.ends_with(".sfasta") {
        format!("{}.sfasta", filename)
    } else {
        filename.to_string()
    }
}

/// Indexes an SFASTA file
pub fn index(filename: &str) -> String {
    // TODO: Run a sanity check on the file first... Make sure it's valid
    // sfasta

    let filesize = metadata(&filename).expect("Unable to open file").len();
    let starting_size = std::cmp::max((filesize / 500) as usize, 1024);

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

    let header: Header = match bincode::deserialize_from(&mut fh) {
        Ok(x) => x,
        Err(_) => panic!("Header missing or malformed in SFASTA file"),
    };

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

    create_index(filename, ids, locations)
}

fn get_index_filename(filename: &str) -> String {
    let filenamepath = Path::new(&filename);
    let filename = Path::new(filenamepath.file_name().unwrap())
        .file_stem()
        .unwrap()
        .to_str()
        .unwrap()
        .to_owned()
        + ".sfai";

    let mut path = filenamepath.parent().unwrap().to_str().unwrap().to_owned();
    if !path.is_empty() {
        path += "/";
    }

    path + &filename
}

fn load_index(filename: &str) -> Option<(Vec<String>, Vec<u64>)> {
    let idx_filename = get_index_filename(filename);
    if !Path::new(&idx_filename).exists() {
        println!("IdxFile does not exist! {} {}", filename, idx_filename);
        return None;
    }

    let (_, _, mut idxfh) = generic_open_file(&idx_filename);
    // let idx: HashMap<String, u64>;
    // let idx: HashMap<String, u64>;
    // idx = bincode::deserialize_from(&mut idxfh).expect("Unable to open Index file");
    let length: u64 = bincode::deserialize_from(&mut idxfh).expect("Unable to read length of index");
    let mut keys: Vec<String> = Vec::with_capacity(length as usize);
    keys = bincode::deserialize_from(&mut idxfh).expect("Unable to read idx keys");
    let mut vals: Vec<u64> = Vec::with_capacity(length as usize);
    vals = bincode::deserialize_from(&mut idxfh).expect("Unable to read idx values");

    Some((keys, vals))
}

fn create_index(filename: &str, ids: Vec<String>, locations: Vec<u64>) -> String {
    let idx: HashMap<String, u64> = ids.into_iter().zip(locations).collect();

    let mut sorted: Vec<_> = idx.into_iter().collect();
    sorted.sort_by(|x,y| x.0.cmp(&y.0));
    let keys: Vec<_> = sorted.iter().map(|x| x.0.clone()).collect();
    let vals: Vec<_> = sorted.iter().map(|x| x.1.clone()).collect();
    
    let output_filename = get_index_filename(filename);

    let out_file = snap::write::FrameEncoder::new(
        File::create(output_filename.clone()).expect("Unable to write to file"),
    );

    let mut out_fh = BufWriter::with_capacity(1024 * 1024, out_file);
    // bincode::serialize_into(&mut out_fh, &idx).expect("Unable to write index");
    
    bincode::serialize_into(&mut out_fh, &(keys.len() as u64)).expect("Unable to write index");
    bincode::serialize_into(&mut out_fh, &keys).expect("Unable to write index");
    bincode::serialize_into(&mut out_fh, &vals).expect("Unable to write index");

    output_filename
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::SeekFrom;

    #[test]
    pub fn convert_fasta_to_sfasta_and_index() {
        let input_filename = "test_data/test_multiple.fna";
        let output_filename = "test_data/test_sfasta_convert_and_index.sfasta";

        convert_fasta_file(input_filename, output_filename);

        let idx_filename = index(output_filename);
        assert!(idx_filename == "test_data/test_sfasta_convert_and_index.sfai");

        load_index(output_filename);
    }

    #[test]
    pub fn test_index() {
        let input_filename = "test_data/test_multiple.fna";
        let output_filename = "test_data/test_index.sfasta";

        convert_fasta_file(input_filename, output_filename);

        let (mut reader, idx) = open_file(output_filename);

        let idx = idx.unwrap();
        let i: Vec<&u64> = idx.1.iter().skip(2).take(1).collect();

        reader
            .seek(SeekFrom::Start(*i[0]))
            .expect("Unable to work with seek API");
        match bincode::deserialize_from::<_, EntryCompressed>(&mut reader) {
            Ok(x) => x,
            Err(_) => panic!("Unable to read indexed SFASTA after jumping"),
        };
    }

    #[test]
    pub fn test_random_mode() {
        let input_filename = "test_data/test_multiple.fna";
        let output_filename = "test_data/test_random.sfasta";

        convert_fasta_file(input_filename, output_filename);

        let mut sequences = Sequences::new(output_filename);
        sequences.set_mode(SeqMode::Random);
        sequences.next();
    }
}
