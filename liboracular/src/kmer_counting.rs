use crossbeam::atomic::AtomicCell;
use once_cell::sync::OnceCell;
use wyhash::wyhash;
use thincollections::thin_vec::ThinVec;
use opinionated::fasta::{complement_nucleotides, capitalize_nucleotides};

use std::sync::{Arc, RwLock};

use std::thread;
use std::thread::Builder;

use std::fs::File;
use std::io::prelude::*;
use std::io::{BufReader, Read, BufRead};
use std::hash::BuildHasherDefault;
use std::collections::HashMap;

use crossbeam::queue::{ArrayQueue, PushError};
use crossbeam::utils::Backoff;

use indicatif::ProgressBar;
use indicatif::ProgressStyle;

use std::time::{Instant};

use serde::{Serialize, Deserialize};

use twox_hash::XxHash;
use fasthash::*;

use std::convert::TryInto;

const STACKSIZE: usize = 256 * 1024 * 1024;  // Stack size (needs to be > BUFSIZE + SEQBUFSIZE)
const WORKERSTACKSIZE: usize = 64 * 1024 * 1024;  // Stack size (needs to be > BUFSIZE + SEQBUFSIZE)

// type Sequences = Vec<Sequence>;
type Sequence = Vec<u8>;

enum ThreadCommand<T> {
    Work(T),
    Terminate,
}

impl ThreadCommand<Sequence> {
    // Consumes the ThreadCommand, which is just fine...
    fn unwrap(self) -> Sequence {
        match self {
            ThreadCommand::Work(x)   => x,
            ThreadCommand::Terminate => panic!("Unable to unwrap terminate command"),
        }
    }
}

pub const MAX_VOCAB: usize = 300_000_000;

pub struct Dict {
    pub wordidx: Vec<AtomicCell<Option<core::num::NonZeroU64>>>,
    pub words:   Vec<OnceCell<ThinVec<u8>>>,
    pub counts:  Vec<AtomicCell<u64>>,
    pub tokens:  AtomicCell<u64>,
    pub size:    AtomicCell<u64>,
    pub entries: AtomicCell<u64>,
}

#[derive(Serialize, Deserialize)]
pub struct FinalDict {
    pub words: HashMap<Vec<u8>, u64, BuildHasherDefault<XxHash>>,
    pub entries: u64,
    pub size: u64,
    pub tokens: u64,
}

impl Dict {

    pub fn convert_to_final(&self) -> FinalDict {

        let mut words: HashMap<Vec<u8>, u64, BuildHasherDefault<XxHash>> = Default::default(); //HashMap::with_capacity(self.size.load() as usize);
        words.reserve(self.size.load() as usize);

        for x in self.words
                    .iter()
                    .filter_map(|x| x.get())
                    .map(|x| x.to_vec()) {
            let id = self.get_id(&x);
            let count = self.counts[self.wordidx[id].load().expect("Error getting wordidx[id]").get() as usize].load();
            words.insert(x, count);
        }

        FinalDict { words, 
                    entries: self.entries.load(), 
                    size: self.size.load(), 
                    tokens: self.tokens.load() 
                }
    }

    pub fn new() -> Dict {

        // Multi-threaded initialization: 24 seconds
        // Single-threaded: 63 seconds...

        /* let mut wordidx = Vec::with_capacity(MAX_VOCAB);
        wordidx.resize_with(MAX_VOCAB, || AtomicCell::new(None));
        // wordidx = (0..MAX_VOCAB).map(|_| AtomicCell::new(None)).collect();
        let mut counts = Vec::with_capacity(MAX_VOCAB);
        // counts = (0..MAX_VOCAB).map(|_| AtomicCell::new(0)).collect();
        counts.resize_with(MAX_VOCAB, || AtomicCell::new(0));
        let mut words = Vec::with_capacity(MAX_VOCAB);
        // words = (0..MAX_VOCAB).map(|_| OnceCell::new()).collect();
        words.resize_with(MAX_VOCAB, OnceCell::new); */


        let wordidx_builder = match Builder::new()
                        .name("WordIdx Builder".into())
                        .spawn(|| 
                        {
                            let mut wordidx = Vec::with_capacity(MAX_VOCAB);
                            wordidx.resize_with(MAX_VOCAB, || AtomicCell::new(None));
                            wordidx
                        })
                    {
                        Ok(x)  => x,
                        Err(y) => panic!("{}", y)
                    };

        let counts_builder = match Builder::new()
                        .name("Counts Builder".into())
                        .spawn(|| 
                        {
                            let mut counts = Vec::with_capacity(MAX_VOCAB);
                            counts.resize_with(MAX_VOCAB, || AtomicCell::new(0));
                            counts
                        })
                    {
                        Ok(x)  => x,
                        Err(y) => panic!("{}", y)
                    };

        let words_builder = match Builder::new()
                        .name("Words Builder".into())
                        .spawn(|| 
                        {
                            let mut words = Vec::with_capacity(MAX_VOCAB);
                            words.resize_with(MAX_VOCAB, OnceCell::new);
                            words
                        })
                    {
                        Ok(x)  => x,
                        Err(y) => panic!("{}", y)
                    };


        let wordidx = wordidx_builder.join().unwrap();
        let counts = counts_builder.join().unwrap();
        let words = words_builder.join().unwrap(); 


        Dict {
            wordidx,
            words,
            counts,
            tokens: AtomicCell::new(0),
            size: AtomicCell::new(1),
            entries: AtomicCell::new(0),
        }
    }

    #[inline]
    fn get_id(&self, kmer: &[u8]) -> usize {
        let rc = self.get_rc(&kmer);

        let mut id = self.calc_hash(&kmer, &rc);
        let mut cur_word = self.words[id].get();
        // let kvec = kmer.to_vec();

        while self.wordidx[id].load() != None 
            && 
            // self.words[id].load() != Some(kmer)
            // &&
            // self.words[id].load() != Some(rc)
            // self.words.read().unwrap()[self.wordidx[id].load().unwrap()].kmer != kmer
            // DashMap below
            //*self.words.index(&id) != kmer
            // &&
            // *self.words.index(&id) != rc
            // self.words.read().unwrap()[self.wordidx[id].load().unwrap()].kmer != rc
            cur_word != None 
            &&
            cur_word.unwrap() != &kmer
            &&
            cur_word.unwrap() != &rc
        {
            id = (id + 1) % MAX_VOCAB;
            cur_word = self.words[id].get();
        }

        id as usize
    }

    #[inline]
    pub fn get_rc(&self, kmer: &[u8]) -> ThinVec<u8> {
        // let mut rc: Vec<u8> = kmer.to_vec();
        let mut rc: ThinVec<u8> = ThinVec::with_capacity(kmer.len());
        rc.extend_from_slice(&kmer);
        complement_nucleotides(&mut rc);
        rc.reverse();
        rc

    }

    #[inline(always)]
    pub fn calc_hash(&self, kmer: &[u8], rc: &[u8]) -> usize {
        // let mut hasher = twox_hash::XxHash64::with_seed(0);
        // hasher.consume(&kmer);
        // hasher.write(&kmer);
        // let x = hasher.finish() as usize % MAX_VOCAB;

        // let mut hasher = twox_hash::XxHash64::with_seed(0);
        // hasher.write(&rc);
        // let y = hasher.finish() as usize % MAX_VOCAB;

        // kmer.iter().map(|x| hasher.write(&x));

        // let rc = self.get_rc(&kmer);

        // SeaHash seems to be faster for this datatype...

        // let x = seahash::hash(&kmer) as usize % MAX_VOCAB;
        // let y = seahash::hash(&rc) as usize % MAX_VOCAB;
        
        // Converting to bits seems slower
        // let val = super::convert_seq_to_bits(&kmer);
        // let  rc = super::bits_rc(val.clone());
        // let x = val.as_slice()[0] as usize % MAX_VOCAB;
        // let y = rc.as_slice()[0] as usize % MAX_VOCAB;

        // Trying FNV hash
        // Fastest so far...
        /* let mut hasher = FnvHasher::with_key(0);
        hasher.write(&kmer);
        let x = hasher.finish() as usize % MAX_VOCAB;

        let mut hasher = FnvHasher::with_key(0);
        hasher.write(&rc);
        let y = hasher.finish() as usize % MAX_VOCAB; */

        // Trying Wy Hash
        // let mut hasher = WyHash::with_seed(3);
        // hasher.write(&kmer);
        // let x = hasher.finish() as usize % MAX_VOCAB;

        // let mut hasher = WyHash::with_seed(3);
        // hasher.write(&rc);
        // let y = hasher.finish() as usize % MAX_VOCAB;

        // Very slow...
        // let x: usize = self.transmute_hash(&kmer) % MAX_VOCAB;
        // let y: usize = self.transmute_hash(&rc) % MAX_VOCAB;

        // Super slow...
        // let x = read_ne_u64(&kmer) as usize % MAX_VOCAB;
        // let y = read_ne_u64(&rc) as usize % MAX_VOCAB;

        let x = wyhash(&kmer, 43_988_123) as usize % MAX_VOCAB;
        let y = wyhash(&rc, 43_988_123) as usize % MAX_VOCAB;

        // Slightly faster than wyhash
        // let x = fasthash::xxh3::hash64(&kmer) as usize % MAX_VOCAB;
        // let y = fasthash::xxh3::hash64(&rc) as usize % MAX_VOCAB;

        // let x = enumerate_sum(&kmer) as usize % MAX_VOCAB;
        // let y = enumerate_sum(&rc) as usize % MAX_VOCAB;

        // Equal to xxh3
        // let x = fasthash::t1ha::hash64(&kmer) as usize % MAX_VOCAB;
        // let y = fasthash::t1ha::hash64(&rc) as usize % MAX_VOCAB;

        std::cmp::min(x,y)
    }

    pub fn add(&self, kmer: &[u8]) {
        self.tokens.fetch_add(1);

        let id = self.get_id(kmer);

        if self.wordidx[id].load() == None {

            let mut word = ThinVec::with_capacity(kmer.len());
            word.extend_from_slice(&kmer);
            if let Err(val) = self.words[id].set(word) {
                // Another thread beat us to it, start over...
                self.add(&val); 
                return; 
            }

            let wordidx = core::num::NonZeroU64::new(self.size.fetch_add(1) as u64);
            
            match self.wordidx[id].compare_and_swap(None, wordidx) {
                None  => (),
                Some(_) => // Collision
                {
                    self.add(kmer);
                    return;
                }
            }; // Points to word

            self.entries.fetch_add(1);
        }
        self.counts[self.wordidx[id].load().unwrap().get() as usize].fetch_add(1);
    }

}

pub fn count_kmers(
    num_threads: usize,
    kmer_size: usize,
    filename: &str) -> Arc<Dict> {

    let now = Instant::now();
    let dict_builder = match Builder::new()
                        .name("Dict Builder".into())
                        .spawn(|| Arc::new(Dict::new()))
                    {
                        Ok(x)  => x,
                        Err(y) => panic!("{}", y)
                    };

    let jobs = Arc::new(AtomicCell::new(0 as usize));
    let seq_queue = Arc::new(ArrayQueue::<ThreadCommand<Sequence>>::new(1024 * 64));
    let rawseq_queue = Arc::new(ArrayQueue::<ThreadCommand<Sequence>>::new(1024));
    let done = Arc::new(RwLock::new(false));
    let generator_done = Arc::new(RwLock::new(false));

    let generator;

    let mut children = Vec::new();

    let filename = filename.to_string();

    for _ in 0..4 {
        let seq_queue = Arc::clone(&seq_queue);
        let rawseq_queue = Arc::clone(&rawseq_queue);
        let jobs = Arc::clone(&jobs);

        let child = match Builder::new()
                        .name("IOWorker".into())
                        .stack_size(WORKERSTACKSIZE)
                        .spawn(move || io_worker_thread(kmer_size, rawseq_queue, seq_queue, jobs)) {
                            Ok(x)  => x,
                            Err(y) => panic!("{}", y)
                        };
        
        children.push(child);
    }

    let dict = dict_builder.join().unwrap();

    {
        let generator_done = Arc::clone(&generator_done);
        let seq_queue = Arc::clone(&seq_queue);
        let rawseq_queue = Arc::clone(&rawseq_queue);
        let jobs = Arc::clone(&jobs);
        let dict = Arc::clone(&dict);
        let mut buffer: Sequence = Vec::with_capacity(1024);

        generator = thread::Builder::new()
                            .name("Generator".to_string())
                            .stack_size(STACKSIZE)
                            .spawn(move||
        {
            let mut seqbuffer: Sequence = Vec::with_capacity(8 * 1024 * 1024); // 8 Mb to start, will likely increase...
            let mut seqlen: usize = 0;

            let file = match File::open(&filename) {
                Err(why) => panic!("Couldn't open {}: {}", filename, why.to_string()),
                Ok(file) => file,
            };

            let pb = ProgressBar::new(file.metadata().unwrap().len());

            let file = BufReader::with_capacity(64 * 1024 * 1024, file);

            pb.set_style(ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:50.cyan/blue} {pos:>7}/{len:7} {eta_precise} eta\n{msg}")
                .progress_chars("█▇▆▅▄▃▂▁  "));

            let fasta: Box<dyn Read> = if filename.ends_with("gz") {
                Box::new(flate2::read::GzDecoder::new(pb.wrap_read(file)))
            } else {
                Box::new(pb.wrap_read(file))
            };

            let mut reader = BufReader::with_capacity(128 * 1024 * 1024, fasta);

            let backoff = Backoff::new();

            while let Ok(bytes_read) = reader.read_until(b'\n', &mut buffer) {
                if bytes_read == 0 {
                    // File is empty, we are done!
                    break;
                }

                match buffer[0] {
                    // 62 is a > meaning we have a new sequence id.
                    62 => {
                        jobs.fetch_add(1 as usize);
                        let wp = ThreadCommand::Work(seqbuffer[..seqlen].to_vec());
                        seqbuffer.clear();
                        seqlen = 0;

                        let mut result = rawseq_queue.push(wp);
                        while let Err(PushError(wp)) = result {
                            result = rawseq_queue.push(wp);
                        }

                        pb.set_message(&format!("{}/1024 {}/65536 {} unique kmers", rawseq_queue.len(), seq_queue.len(), dict.entries.load()));
                    },

                    // Anything else is likely sequence we need...
                    _  => {
                        //seqbuffer.resize(seqlen + bytes_read - 1, 0);
                        // seqbuffer[seqlen..seqlen + bytes_read - 1].copy_from_slice(&buffer[0..bytes_read - 1]);
                        let slice_end = bytes_read.saturating_sub(1);
                        seqbuffer.extend_from_slice(&buffer[0..slice_end]);
                        seqlen = seqlen.saturating_add(slice_end);
                        // println!("{:#?}", String::from_utf8(seqbuffer[0..seqlen].to_vec()).unwrap());
                    }
                }

            buffer.clear();
            }

            *generator_done.write().unwrap() = true;
        }).unwrap();
    }

    for _ in 0..num_threads {
        let dict = Arc::clone(&dict);
        let seq_queue = Arc::clone(&seq_queue);
        let jobs = Arc::clone(&jobs);

        let child = match Builder::new()
                        .name("Worker".into())
                        .stack_size(WORKERSTACKSIZE)
                        .spawn(move || kmer_counter_worker_thread(kmer_size, dict, seq_queue, jobs)) {
                            Ok(x)  => x,
                            Err(y) => panic!("{}", y)
                        };
        
        children.push(child);
    }

    let backoff = Backoff::new();
    while !*generator_done.read().unwrap() {
        backoff.snooze();
    }

    generator.join().expect("Unable to join generator thread...");

    // IO Generator is done, shut down the IO workers...
    for _ in 0..4 {
        match rawseq_queue.push(ThreadCommand::Terminate) {
            Ok(_) => (),
            Err(x) => panic!("Unable to send command... {:#?}", x)
        }
    }

    println!("Waiting for seq_queue to be empty...Currently at: {}", seq_queue.len());

    while !seq_queue.is_empty() {
        backoff.snooze();
    }

    println!("Seq queue is empty, sending terminate command...");

    while jobs.load() > 0 {
        backoff.snooze();
        // println!("{}", jobs.load());
    }

    for _ in 0..num_threads {
        match seq_queue.push(ThreadCommand::Terminate) {
            Ok(_) => (),
            Err(x) => panic!("Unable to send command... {:#?}", x)
        }
    }

    println!("Terminate commands sent, joining worker threads");

    for child in children {
        match child.join() {
            Ok(_) => (),
            Err(x) => panic!("Error joining worker thread... {:#?}", x)
        }
    }

    println!("Worker threads joined, getting dictionary stats...");

    *done.write().unwrap() = true;

    dict

}

fn io_worker_thread(
    kmer_size: usize,
    rawseq_queue: Arc<ArrayQueue<ThreadCommand<Sequence>>>,
    seq_queue: Arc<ArrayQueue<ThreadCommand<Sequence>>>,
    jobs: Arc<AtomicCell<usize>>) {
    
    let backoff = Backoff::new();

    loop {
        if let Ok(command) = rawseq_queue.pop() {
            // We got the work packet!

            // We are finished, end the thread...
            if let ThreadCommand::Terminate = command {
                break;
            }

            let mut rawseq = command.unwrap();
            jobs.fetch_sub(1);

            capitalize_nucleotides(&mut rawseq);
            let coords = super::utils::get_good_sequence_coords(&rawseq);
            
            for (start_coords, end_coords) in coords {

                if (end_coords - start_coords) >= kmer_size {

                    jobs.fetch_add(1 as usize);
                    let wp = ThreadCommand::Work(rawseq[start_coords..end_coords].to_vec());

                    let mut result = seq_queue.push(wp);
                    while let Err(PushError(wp)) = result {
                        backoff.snooze();
                        result = seq_queue.push(wp);
                    }

                }
                
            }

        } else {
            backoff.snooze();
            backoff.reset();
        }
    }


}

fn kmer_counter_worker_thread (
    kmer_size: usize, 
    dict: Arc<Dict>, 
    seq_queue: Arc<ArrayQueue<ThreadCommand<Sequence>>>,
    jobs: Arc<AtomicCell<usize>>) {

    let backoff = Backoff::new();

    loop {
        /* while seq_queue.is_empty() {
            backoff.snooze();
        } */

        // Even when not empty, not guaranteed to be the thread to grab the work first...
        if let Ok(command) = seq_queue.pop() {
            // We got the work packet!

            // We are finished, end the thread...
            if let ThreadCommand::Terminate = command {
                break;
            }          

            let rawseq = command.unwrap();
            jobs.fetch_sub(1);

            for i in 0..kmer_size {
                rawseq[i..].chunks_exact(kmer_size).for_each(|x| { 
                    // if bytecount::count(&x, b'N') < 3 {
                        dict.add(&x);
                    // }
                });
            }

        } else {
            backoff.snooze();
            backoff.reset();
        }

    }
}

#[cfg(test)]
mod test {
    use crossbeam::atomic::AtomicCell;
    use once_cell::sync::OnceCell;
    use thincollections::thin_vec::ThinVec;
    use std::mem;

    #[test]
    fn is_NOT_lock_free_NonZeroU32() {
        assert_eq!(AtomicCell::<Option<core::num::NonZeroU32>>::is_lock_free(), false);
    }

    #[test]
    fn is_NOT_lock_free_u32() {
        assert_eq!(AtomicCell::<u32>::is_lock_free(), false);
    }

    #[test]
    fn is_lock_free_usize() {
        assert_eq!(AtomicCell::<usize>::is_lock_free(), true);
    }

    #[test]
    fn is_lock_free_u64() {
        assert_eq!(AtomicCell::<u64>::is_lock_free(), true);
    }

    #[test]
    fn is_lock_free_i64() {
        assert_eq!(AtomicCell::<i64>::is_lock_free(), true);
    }

    #[test]
    fn size_of() {
        assert_eq!(mem::size_of::<OnceCell<ThinVec<u8>>>(), 16);
    }

    #[test]
    fn is_lock_free_NonZeroU64() {
        assert_eq!(AtomicCell::<Option<core::num::NonZeroU64>>::is_lock_free(), true);
    }

    


}