use crossbeam::atomic::AtomicCell;
use once_cell::sync::OnceCell;
use wyhash::wyhash;
use thincollections::thin_vec::ThinVec;
use opinionated::fasta::{complement_nucleotides, capitalize_nucleotides};

use std::sync::{Arc};

use std::thread::Builder;

use std::hash::BuildHasherDefault;
use std::collections::HashMap;

use crossbeam::queue::{ArrayQueue, PushError};
use crossbeam::utils::Backoff;

use serde::{Serialize, Deserialize};

use twox_hash::XxHash;

use crate::threads::{sequence_generator, Sequence, ThreadCommand};

const STACKSIZE: usize = 256 * 1024 * 1024;  // Stack size (needs to be > BUFSIZE + SEQBUFSIZE)
const WORKERSTACKSIZE: usize = 64 * 1024 * 1024;  // Stack size (needs to be > BUFSIZE + SEQBUFSIZE)

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

        // Overkill, but multi-threading gives it a slight speed boost
        // Higher in dev...

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

        let x = wyhash(&kmer, 43_988_123) as usize % MAX_VOCAB;
        let y = wyhash(&rc, 43_988_123) as usize % MAX_VOCAB;

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
                Some(_) => // Collision... unlikely to happen, but it could...
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

    let dict_builder = match Builder::new()
                        .name("Dict Builder".into())
                        .spawn(|| Arc::new(Dict::new()))
                    {
                        Ok(x)  => x,
                        Err(y) => panic!("{}", y)
                    };

    let (seq_queue, jobs, generator_done, generator, mut children) 
        = sequence_generator(kmer_size, &filename);

    let dict = dict_builder.join().unwrap();

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

    dict

}

fn kmer_counter_worker_thread (
    kmer_size: usize, 
    dict: Arc<Dict>, 
    seq_queue: Arc<ArrayQueue<ThreadCommand<Sequence>>>,
    jobs: Arc<AtomicCell<usize>>) {

    let backoff = Backoff::new();

    loop {
        if let Ok(command) = seq_queue.pop() {

            // We are finished, end the thread...
            if let ThreadCommand::Terminate = command {
                return;
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
    fn is_not_lock_free_non_zero_u32() {
        assert_eq!(AtomicCell::<Option<core::num::NonZeroU32>>::is_lock_free(), false);
    }

    #[test]
    fn is_not_lock_free_u32() {
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
    fn is_lock_free_non_zero_u64() {
        assert_eq!(AtomicCell::<Option<core::num::NonZeroU64>>::is_lock_free(), true);
    }

    


}