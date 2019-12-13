use crossbeam::atomic::AtomicCell;
use once_cell::sync::OnceCell;
use thincollections::thin_vec::ThinVec;
use opinionated::fasta::{complement_nucleotides};

use std::sync::{Arc};

use std::mem;

use std::thread::Builder;

use std::hash::BuildHasherDefault;
use std::collections::HashMap;

use crossbeam::queue::{ArrayQueue};
use crossbeam::utils::Backoff;

use serde::{Serialize, Deserialize};

use t1ha::{t1ha0};

use fnv::FnvHashMap;
use fnv::FnvHasher;

use indicatif::ProgressBar;
use indicatif::ProgressStyle;

use crate::threads::{sequence_generator, Sequence, ThreadCommand};
use crate::kmer_hasher::{kmerhash, kmerhash_smallest, calc_rc, hash4, convert_to_kmer};

// const STACKSIZE: usize = 256 * 1024 * 1024;  // Stack size (needs to be > BUFSIZE + SEQBUFSIZE)
const WORKERSTACKSIZE: usize = 16 * 1024 * 1024;  // Stack size (needs to be > BUFSIZE + SEQBUFSIZE)

/* pub const BUCKET_SIZE: usize = 31;
pub const MAX_VOCAB:   usize = 1 << BUCKET_SIZE; // 2 147 483 648
pub const HALF_VOCAB:  usize = 1 << (BUCKET_SIZE - 1);
pub const MAX_MASK:    usize = MAX_VOCAB - 1;
pub const HALF_MASK:   usize = HALF_VOCAB - 1; */

// TODO: Switch words to u64, to make use of faster "hash" and faster RC computation
pub struct Dict {
    pub wordidx: Vec<AtomicCell<Option<core::num::NonZeroU64>>>,
    pub words:   Vec<OnceCell<u64>>,
    pub counts:  Vec<AtomicCell<u64>>,
    pub tokens:  AtomicCell<u64>,
    pub size:    AtomicCell<u64>,
    pub entries: AtomicCell<u64>,
    pub k: u64,
    pub bucket_size: usize,
    pub max_vocab: usize,
    pub half_vocab: usize,
    pub max_mask: usize,
    pub half_mask: usize
}

#[derive(Serialize, Deserialize)]
pub struct FinalDict {
    // pub words: HashMap<Vec<u8>, u64, BuildHasherDefault<XxHash>>,
    pub words: HashMap<String, u64, BuildHasherDefault<FnvHasher>>,
    pub entries: u64,
    // pub size: u64,
    pub tokens: u64,
    pub k: u64,
}

impl Dict {

    // TODO: Make multi-threaded, this is the slow part now...
    pub fn convert_to_final(&self) -> FinalDict {
        
        let mut words: HashMap<String, u64, _> = FnvHashMap::default();
        words.reserve(self.size.load() as usize);

        let mut tokens = 0;

        /* for x in self.words
                    .iter()
                    .filter_map(|x| x.get()) { */
            // let id = self.get_id(*x);
        for (id, countid) in self.wordidx
                    .iter()
                    .enumerate()
                    .filter(|(_, countid_)| countid_.load() != None)
                    .filter(|(_, countid_)| countid_.load().unwrap().get() != 0)
                    .map(|(id, countid_)| (id, countid_.load().unwrap().get() as usize))
                    {
                    //.filter_map(|x| x.load()) {
                    //.map(|x| x.get() as usize) {

            let x = match self.words[id].get() {
                Some(x) => x,
                None    => continue //panic!("Not found: {} {}", id, self.counts[id].load())
            };
            // let count = self.counts[self.wordidx[id].load().expect("Error getting wordidx[id]").get() as usize].load();
            let count = self.counts[countid].load();

            // let xstr = unsafe { String::from_utf8_unchecked(x) };
            // let rcstr = unsafe { String::from_utf8_unchecked(rc.to_vec()) };

            let xstr = convert_to_kmer(self.k as usize, *x);
            let rcstr = convert_to_kmer(self.k as usize, calc_rc(self.k as usize, *x));

            *words.entry(xstr).or_insert(0) += count;
            *words.entry(rcstr).or_insert(0) += count;
            tokens += count + count;
        }

        FinalDict { 
                    entries: words.keys().len() as u64, 
                    words, 
                    tokens: tokens,
                    k: self.k,
                }
    }

    pub fn new(k: usize, bucket_size: usize) -> Dict {

        assert!(bucket_size < 64, "Maximum bucket size must be <64");

        let max_vocab = 1 << bucket_size;
        let half_vocab = 1 << (bucket_size - 1);
        let max_mask = max_vocab - 1;
        let half_mask = half_vocab - 1;

        let (wordidx_builder, counts_builder, words_builder);

        // Overkill, but multi-threading gives it a slight speed boost
        // Higher in dev...

        {
            let max_vocab = max_vocab.clone();
            wordidx_builder = match Builder::new()
                        .name("WordIdx Builder".into())
                        .spawn(move || 
                        {
                            let mut wordidx = Vec::with_capacity(max_vocab);
                            wordidx.resize_with(max_vocab, || AtomicCell::new(None));
                            wordidx
                        })
                    {
                        Ok(x)  => x,
                        Err(y) => panic!("{}", y)
                    };
        }

        {
            let max_vocab = max_vocab.clone();
            counts_builder = match Builder::new()
                        .name("Counts Builder".into())
                        .spawn(move || 
                        {
                            let mut counts = Vec::with_capacity(max_vocab);
                            counts.resize_with(max_vocab, || AtomicCell::new(0));
                            counts
                        })
                    {
                        Ok(x)  => x,
                        Err(y) => panic!("{}", y)
                    };
        }

        {
            let max_vocab = max_vocab.clone();
    
            words_builder = match Builder::new()
                        .name("Words Builder".into())
                        .spawn(move || 
                        {
                            let mut words = Vec::with_capacity(max_vocab);
                            words.resize_with(max_vocab, OnceCell::new);
                            words
                        })
                    {
                        Ok(x)  => x,
                        Err(y) => panic!("{}", y)
                    };

        }


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
            k: k as u64,
            bucket_size,
            max_vocab,
            half_mask,
            half_vocab,
            max_mask
        }

/* pub const BUCKET_SIZE: usize = 31;
pub const MAX_VOCAB:   usize = 1 << BUCKET_SIZE; // 2 147 483 648
pub const HALF_VOCAB:  usize = 1 << (BUCKET_SIZE - 1);
pub const MAX_MASK:    usize = MAX_VOCAB - 1;
pub const HALF_MASK:   usize = HALF_VOCAB - 1; */

    }

    #[inline]
    fn get_id(&self, kmer: u64) -> usize {
        // let rc = self.get_rc(&kmer);

        // let hash = kmerhash_smallest(&kmer) as usize;
        //let mut id = self.calc_hash(&kmer); // , &rc
        let kmer_as_u8: [u8; 8] = unsafe { mem::transmute(kmer) };

        // This helps distribute the hashes even better, speeding it up!
        let hash = t1ha0(&kmer_as_u8, self.max_vocab as u64) as usize;

        // let mut id = hash % HALF_VOCAB;
        let mut id = (hash ^ self.half_mask) & self.half_mask;
        let mut cur_word = self.words[id].get();
        let mut retry: usize = 0;

        while self.wordidx[id].load() != None 
            && 
            cur_word != None 
            &&
            cur_word.unwrap() != &kmer
            // &&
            // cur_word.unwrap() != &rc
        {
            retry = retry.saturating_add(1);
            // id = (hash.wrapping_add(retry)) % MAX_VOCAB;
            id = (hash.wrapping_add(retry) ^ self.max_mask) & self.max_mask;
            if retry > 5_000_000 {
              assert!(self.entries.load() < self.max_vocab as u64, "More entries than MAX_VOCAB!");
              // println!("Error: More than 100,000 tries... Tokens: {} Entries: {}", self.tokens.load(), self.entries.load());
              // println!("Retry: {} Kmer is: {} ID is: {}", retry, String::from_utf8(kmer.to_vec()).unwrap(), id);
            }
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
    pub fn calc_hash(&self, kmer: &[u8]) -> usize { // , rc: &[u8]

        // let x = wyhash(&kmer, 43_988_123) as usize % MAX_VOCAB;
        // let y = wyhash(&rc, 43_988_123) as usize % MAX_VOCAB;

        //let x = t1ha0(&kmer, 42_988_123) as usize % MAX_VOCAB;
        //let y = t1ha0(&rc, 42_988_123) as usize % MAX_VOCAB;

        //std::cmp::min(x,y)
        kmerhash_smallest(kmer) as usize % self.max_vocab
    }

    /*
     *  wordidx[id] = idsize
     *    words[id] = kmer
     * 
     *  counts[idsize]
     * 
     */

    pub fn add(&self, kmer: u64) {
        self.tokens.fetch_add(1);

        let id = self.get_id(kmer);
        // let id = kmer as usize;

        if self.wordidx[id].load() == None {
            let word = kmer;
            // let mut word = ThinVec::with_capacity(kmer.len());
            // word.extend_from_slice(&kmer);
            if let Err(val) = self.words[id].set(word) {
                // Another thread beat us to it, start over...
                self.tokens.fetch_sub(1);
                self.add(val); 
                return;
            }

            let wordidx = core::num::NonZeroU64::new(self.size.fetch_add(1) as u64);

            match self.wordidx[id].compare_and_swap(None, wordidx) {
                None  => (),
                Some(_) => // Collision... unlikely to happen, but it could...
                {
                    self.tokens.fetch_sub(1);
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

    let k = kmer_size.clone();
    let dict_builder = match Builder::new()
                        .name("Dict Builder".into())
                        .spawn(move || Arc::new(Dict::new(k, 29))) // TODO: Make argument-- 31 for nt, 28 for Vvulg
                    {
                        Ok(x)  => x,
                        Err(y) => panic!("{}", y)
                    };

    let (seq_queue, jobs, generator_done, generator, mut children) 
        = sequence_generator(kmer_size, &filename, Arc::new("".to_string()), 1_u64);

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
    backoff.snooze();
    backoff.snooze();
    backoff.snooze();

    let pb = ProgressBar::new_spinner();
    pb.enable_steady_tick(120);
    pb.set_style(
        ProgressStyle::default_spinner()
            // For more spinners check out the cli-spinners project:
            // https://github.com/sindresorhus/cli-spinners/blob/master/spinners.json
            .tick_strings(&[
                "▹▹▹▹▹",
                "▸▹▹▹▹",
                "▹▸▹▹▹",
                "▹▹▸▹▹",
                "▹▹▹▸▹",
                "▹▹▹▹▸",
                "▪▪▪▪▪",
            ])
            .template("{spinner:.blue} {msg}"),
    );

    while !*generator_done.read().unwrap() {
        backoff.snooze();
    }

    pb.set_message("Joining generator");
    generator.join().expect("Unable to join generator thread...");
    pb.set_message("Generator joined");

    while !seq_queue.is_empty() {
        pb.set_message(&format!("{} {} {} Waiting for Seq Queue to be empty", dict.size.load(), dict.tokens.load(), jobs.load()));
        backoff.snooze();
        backoff.snooze();
        backoff.snooze();
        backoff.snooze();
        backoff.snooze();
        backoff.snooze();
    }

    while jobs.load() > 0 {
        backoff.snooze();
        pb.set_message(&format!("{} {} {} Finishing Final Jobs", dict.size.load(), dict.tokens.load(), jobs.load()));
        backoff.snooze();
        backoff.snooze();
        backoff.snooze();
        backoff.snooze();
        backoff.snooze();
        backoff.snooze();

        // println!("{}", jobs.load());
    }

    for _ in 0..num_threads {
        pb.set_message(&format!("{} {} {} Finishing Final Jobs", dict.size.load(), dict.tokens.load(), jobs.load()));
        match seq_queue.push(ThreadCommand::Terminate) {
            Ok(_) => (),
            Err(x) => panic!("Unable to send command... {:#?}", x)
        }
    }

    pb.set_message(&format!("{} {} {} Finishing Final Jobs", dict.size.load(), dict.tokens.load(), jobs.load()));

    for child in children {
        pb.set_message(&format!("{} {} {} Finishing Final Jobs", dict.size.load(), dict.tokens.load(), jobs.load()));
        match child.join() {
            Ok(_) => (),
            Err(x) => panic!("Error joining worker thread... {:#?}", x)
        }
    }

    pb.set_message(&format!("Finished counting kmers"));
    pb.finish();
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

            // TODO: Disable for nt... probably...
            //for i in 0..kmer_size {
            // for i in 0..kmer_size {
                // let i = 0;

                // 861.82user 167.09system 6:27.64elapsed 265%CPU (0avgtext+0avgdata 54289764maxresident)k
                // k=21, Vvulg.fna, t=16
                rawseq
                    .chunks_exact(kmer_size)
                    .for_each(|chunk| dict.add(kmerhash_smallest(chunk)));

                // TODO: Handle remainder (could still be > kmer_size)

                // 822.36user 180.69system 7:04.79elapsed 236%CPU (0avgtext+0avgdata 54289232maxresident)k
                // k=21, Vvulg.fna, t=16
                /* for chunk in rawseq.chunks_exact(4 * kmer_size) {
                    // TODO: Probably don't need to iterate this.... just do direct slices...
                    let kmers = chunk.chunks_exact(kmer_size).collect::<Vec<&[u8]>>();
                    // dict.add(kmerhash_smallest(kmers[0]));
                    // dict.add(kmerhash_smallest(kmers[1]));
                    // dict.add(kmerhash_smallest(kmers[2]));
                    // dict.add(kmerhash_smallest(kmers[3]));
                    let hashes = hash4(kmers[0], kmers[1], kmers[2], kmers[3]);
                    dict.add(hashes.0);
                    dict.add(hashes.1);
                    dict.add(hashes.2);
                    dict.add(hashes.3);
                } */
            // }
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

    use crate::kmer_counting::{*};
    use crate::kmer_hasher::{hash4};

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

    #[test]
    fn dictionary_generator() {
        let k = 21;
        let dict_builder = match Builder::new()
                        .name("Dict Builder".into())
                        .spawn(move || Arc::new(Dict::new(k, 8)))
                    {
                        Ok(x)  => x,
                        Err(y) => panic!("{}", y)
                    };
        let dict = dict_builder.join().unwrap();

        let k = (b"ACTCACGATCACGATACAAAN", 
                 b"TCAGTCACTAGCATACAACTC",
                 b"ACGATCACGATACAAANNNNN",
                 b"TGACTANCATCANTACTTGGT");

        let krc = (b"NTTTGTATCGTGATCGTGAGT",
                 b"GAGTTGTATGCTAGTGACTGA",
                 b"NNNNNTTTGTATCGTGATCGT",
                 b"ACCAAGTANTGATGNTAGTCA");

        let hashes   = hash4(k.0, k.1, k.2, k.3);
        let hashesrc = hash4(krc.0, krc.1, krc.2, krc.3);

        dict.add(hashes.0);
        assert_eq!(dict.tokens.load(), 1, "Wrong number of tokens after first add");
        dict.add(hashes.0);
        dict.add(hashes.0);
        dict.add(hashes.0);
        dict.add(hashes.0);
        dict.add(hashes.0);
        assert_eq!(dict.tokens.load(), 6, "Wrong number of tokens after first add");
        assert_eq!(dict.entries.load(), 1, "Invalid number of entries");
        dict.add(hashesrc.0);
        assert_eq!(dict.tokens.load(), 7, "Wrong number of tokens after first add");
        assert_eq!(dict.entries.load(), 1, "Invalid number of entries");

        dict.add(hashes.1);
        assert_eq!(dict.tokens.load(), 8, "Wrong number of tokens after first add");
        assert_eq!(dict.entries.load(), 2, "Invalid number of entries");
    }


    


}
