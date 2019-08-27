use seahash;
use crossbeam::atomic::AtomicCell;
use crossbeam::utils::Backoff;
use std::sync::{Arc, RwLock, RwLockWriteGuard};

// use std::sync::RwLock;

                             // 224 757 086 nt unique kmers at k=13
pub const MAX_VOCAB: usize = 50000000; //real
// pub const MAX_VOCAB: usize = 20000000; // Debugging...

pub struct Dict {
    wordidx: RwLock<Vec<Option<usize>>>, // Randomish, based on hash
    rcidx:   RwLock<Vec<Option<usize>>>, // Given a standard word ID, get the ID of the reverse complement
    words:   RwLock<Vec<Entry>>, // Sequential
    pub tokens:  AtomicCell<usize>,
    pub size:    AtomicCell<u32>,
    pub entries: AtomicCell<usize>,
    minimum_threshold: AtomicCell<u32>,
    maximum_threshold: AtomicCell<u32>,
    ignore_table: RwLock<Vec<bool>>,
    pct75: u32,
    discard_table: RwLock<Vec<usize>>
}

#[derive(Debug)]
pub struct Entry {
    pub kmer: Vec<u8>,
    pub count: AtomicCell<u32>,
}

impl Dict {
    pub fn new() -> Dict {
        let wordidx = RwLock::new(Vec::new());
        wordidx.write().unwrap().resize_with(MAX_VOCAB, || None);

        let rcidx = RwLock::new(Vec::new());
        rcidx.write().unwrap().resize_with(MAX_VOCAB, || None);

        Dict {
            // vec![None; MAX_VOCAB],
            wordidx: wordidx,
            rcidx: rcidx,
            words: RwLock::new(Vec::with_capacity(MAX_VOCAB / 2)),
            ignore_table: RwLock::new(vec![false; MAX_VOCAB]),
            tokens: AtomicCell::new(0),
            size: AtomicCell::new(0),
            entries: AtomicCell::new(0),
            minimum_threshold: AtomicCell::new(1),
            maximum_threshold: AtomicCell::new(std::u16::MAX as u32),
            pct75: (0.75 * MAX_VOCAB as f64) as u32,
            discard_table: RwLock::new(Vec::new())
        }
    }

    #[inline]
    fn get_id(&self, hash: usize, kmer: &[u8]) -> usize {
        let mut id = hash;
        let wordidx_r = self.wordidx.read().unwrap();
        let words_r = self.words.read().unwrap();

        while wordidx_r[id] != None 
            && 
            words_r[wordidx_r[id].unwrap()].kmer != kmer {
            id = (id + 1) % MAX_VOCAB;
        }

        id as usize
    }

    #[inline]
    fn get_id_with_lock(&self, 
        lock: &RwLockWriteGuard<Vec<Option<usize>>>, 
        words_lock: &RwLockWriteGuard<Vec<Entry>>,
        hash: usize, kmer: &[u8]) 
            -> usize {

        let mut id = hash;
        while lock[id] != None 
            && 
            words_lock[lock[id].unwrap()].kmer != kmer {
            id = (id + 1) % MAX_VOCAB;
        }

        id as usize
    }


    /* fn get_id_prev(&self, kmer: &[u8]) -> usize {
        let mut id = (seahash::hash(kmer) as usize) % MAX_VOCAB;
        while self.wordidx[id] != None && self.words[self.wordidx[id].unwrap()].kmer != kmer {
            id = (id + 1) % MAX_VOCAB;
        }

        id as usize
    } */

    //pub fn add(&mut self, kmer: &[u8]) -> usize {
    pub fn add(&self, hash: usize, kmer: Vec<u8>) {
        self.tokens.fetch_add(2);

        if *self.ignore_table.read().unwrap().get(hash).unwrap() {
            return ()
        }

        let id = self.get_id(hash, &kmer);

        if self.wordidx.read().unwrap()[id] == None {
            // self.entrieself.wordidxs += 2;
            self.entries.fetch_add(2);

            let mut wordidx_w = self.wordidx.write().unwrap();
            let mut rcidx_w = self.rcidx.write().unwrap();
            let mut words_w = self.words.write().unwrap();
            
            let wordidx = match self.discard_table.write().unwrap().pop() {
                None => {
                    let id = Some(self.size.load() as usize);
                    self.size.fetch_add(1);
                    words_w.push(Entry { kmer: kmer.to_vec(), count: AtomicCell::new(1), });
                    id
                },
                Some(x) => {
                    words_w[x] = Entry { kmer: kmer.to_vec(), count: AtomicCell::new(1), };
                    Some(x)
                },
            };

            wordidx_w[id] = wordidx; // Points to word

            // Assuming if the word doesn't exist then it's RC doesn't exist...
            let mut rc = kmer.to_vec();
            super::complement_nucleotides(&mut rc);
            rc.reverse();
            let rcid = self.get_id_with_lock(&wordidx_w, &words_w, seahash::hash(&rc) as usize % MAX_VOCAB, &rc);

            rcidx_w[rcid] = wordidx;

            // if wordidx_w[rcid] == None {
                
                let rcwordidx = match self.discard_table.write().unwrap().pop() {
                    None => {
                        let id = Some(self.size.load() as usize);
                        self.size.fetch_add(1);
                        words_w.push(Entry { kmer: rc, count: AtomicCell::new(1), });
                        id
                    }
                    Some(x) => {
                        words_w[x] = Entry { kmer: rc, count: AtomicCell::new(1), };
                        Some(x)
                    },
                };

                wordidx_w[rcid] = rcwordidx;
                rcidx_w[id] = rcwordidx; // Points to rc idx of word

            // }

            drop(wordidx_w);
            drop(rcidx_w);
            drop(words_w);

            // RC's are either +1 or -1 so we can't simply add one...

            // k of 13 gives LOTS of singletons...
            if self.size.load() >= self.pct75 {
                println!("Pruning...which probably doesn't work anymore...");
                // self.minimum_threshold += 1;
                self.minimum_threshold.fetch_add(1);
                
                let mut singletons: Vec<_> = self.words.read().unwrap()
                    .iter()
                    .enumerate()
                    .filter(|(_ , x)| x.count.load() < self.minimum_threshold.load())
                    .map(|(y, _)| y)
                    .collect();

                println!("");
                println!("We are at 75%!!!!");
                println!("We are at 75%!!!! Current total: {} {}", self.words.read().unwrap().len(), singletons.len());
                println!("Highest count is: {}", self.words.read().unwrap().iter().map(|x| x.count.load()).max().unwrap());
                println!("Potentially informative: {}", self.words.read().unwrap().len() - singletons.len());
                println!("");

                singletons.sort_unstable();
                singletons.reverse();
                let removed = singletons.len();

                for k in singletons {
                    let v = &self.words.read().unwrap()[k];
/*                    let id = self.get_id(
                        (seahash::hash(&v.kmer) as usize % MAX_VOCAB) as usize, 
                        &v.kmer); */
                    self.ignore_table.write().unwrap()[(seahash::hash(&v.kmer) as usize % MAX_VOCAB) as usize] = true;
                    self.discard_table.write().unwrap().push(k);
                }

                self.entries.fetch_sub(removed);

                let toomany: Vec<_> = self.words.read().unwrap()
                    .iter()
                    .enumerate()
                    .filter(|(_ , x)| x.count.load() > self.maximum_threshold.load())
                    .map(|(y, _)| y)
                    .collect();

                for k in toomany {
                    let v = &self.words.read().unwrap()[k];
/*                    let id = self.get_id(
                        (seahash::hash(&v.kmer) as usize % MAX_VOCAB) as usize, 
                        &v.kmer); */
                    self.ignore_table.write().unwrap()[(seahash::hash(&v.kmer) as usize % MAX_VOCAB) as usize] = true;
                    self.discard_table.write().unwrap().push(k);
                }

                self.entries.fetch_sub(removed);

                println!("");
                println!("Done with the pruning...");
                println!("{}", self.size.load());
                println!("");
            }

        } else {
            // println!("{}", String::from_utf8(kmer.to_vec()).unwrap());
            // self.words[self.wordidx[id].unwrap()].count += 1;
            // self.words[self.rcidx[id].unwrap()].count += 1;

            // let mut cnt = self.words.read().unwrap()[self.wordidx.read().unwrap()[id].unwrap()].read().unwrap().count.saturating_add(1);
            // self.words.write().unwrap()[self.wordidx.read().unwrap()[id].unwrap()].write().unwrap().count = cnt;

            let mut added = false;

            let backoff = Backoff::new();
            
            while !added {

                let words_r = self.words.try_read();
                let wordidx_r = self.wordidx.try_read();
                let rcidx_r = self.rcidx.try_read();

                match (words_r, wordidx_r, rcidx_r) {
                    (Ok(words_r), Ok(wordidx_r), Ok(rcidx_r)) => {
                        words_r[wordidx_r[id].unwrap()].count.fetch_add(1);
                        words_r[rcidx_r[id].unwrap()].count.fetch_add(1);
                        added = true;
                    },
                    _ => {}
                }

                if added {
                    break;
                }

                backoff.spin();
            }
           

/*            match self.rcidx.read().unwrap()[id] {
                None => {
                    println!("{:#?}", self.rcidx.read().unwrap()[id]);
                    println!("Words Len: {}", self.words.read().unwrap().len());
                    println!("{:#?}", self.words.read().unwrap()[self.wordidx.read().unwrap()[id].unwrap()].read().unwrap());
                    println!("{}", self.wordidx.read().unwrap()[id].unwrap());
                    panic!("Is none: self.rcidx[id].load(): {}", id);
                },
                Some(_) => ()
            }; */

/*            match self.words.read().unwrap()[self.rcidx.read().unwrap()[id].unwrap()].read() {
                Err(y) => panic!("Error 2: {} {}", id, y),
                Ok(_) => ()
            }; */

            // println!("Word entry: {:#?}", self.words.read().unwrap()[self.rcidx[id].load().unwrap()]);

            // cnt = self.words.read().unwrap()[self.rcidx.read().unwrap()[id].unwrap()].read().unwrap().count.saturating_add(1);
            // self.words.write().unwrap()[self.rcidx.read().unwrap()[id].unwrap()].write().unwrap().count = cnt;
            // self.words.read().unwrap()[self.rcidx.read().unwrap()[id].unwrap()].read().unwrap().count.fetch_add(1);
        }

        // Regular addition... (+ symbol)
        // 80.90user 22.52system 1:11.13elapsed 145%CPU (0avgtext+0avgdata 12235840maxresident)k

        // wrapping_add
        // 75.55user 20.42system 1:04.16elapsed 149%CPU (0avgtext+0avgdata 12240500maxresident)k

        // saturating_add
        // 75.55user 20.42system 1:04.16elapsed 149%CPU (0avgtext+0avgdata 12240500maxresident)k

        // self.wordidx[id].unwrap()
    }

/*    pub fn add_previous(&mut self, kmer: &[u8]) {
        self.tokens += 2;
        let id = self.get_id(&kmer);
        if self.wordidx[id] == None {
            self.words.push(Entry { kmer: kmer.to_vec(), count: 1, });
            let wordidx = Some(self.size as usize);
            self.size += 1;
            self.wordidx[id] = wordidx; // Points to word

            // Assuming if the word doesn't exist then it's RC doesn't exist...
            let mut rc = kmer.to_vec();
            super::complement_nucleotides(&mut rc);
            rc.reverse();
            let rcid = self.get_id(&rc);

            self.rcidx[rcid] = wordidx;

            if self.wordidx[rcid] == None {
                self.words.push(Entry { kmer: rc, count: 1, });
                let rcwordidx = Some(self.size as usize);
                self.size += 1; 

                self.wordidx[rcid] = rcwordidx;
                self.rcidx[id] = rcwordidx; // Points to rc idx of word
            }            
            // RC's are either +1 or -1 so we can't simply add one...

            if (self.size % 10000000) == 0 {
                println!("{}", self.size);
            }

        } else {
            // println!("{}", String::from_utf8(kmer.to_vec()).unwrap());
            self.words[self.wordidx[id].unwrap()].count += 1;
            self.words[self.rcidx[id].unwrap()].count += 1;
        }
        // self.wordidx[id].unwrap()
    } */

    /*
    pub fn add_to_id(&mut self, id: usize) {
        self.words[id].count += 1;
    } */

}