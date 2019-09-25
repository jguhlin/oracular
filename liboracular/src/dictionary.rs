use crossbeam::atomic::AtomicCell;
use once_cell::sync::OnceCell;
use wyhash::wyhash;
use thincollections::thin_vec::ThinVec;

// use num_traits::int::{u64, u8, usize};

// use std::sync::RwLock;
// TODO: Try out https://github.com/matklad/once_cell/ instead of DashMap
                             // 224 757 086 nt unique kmers at k=13
pub const MAX_VOCAB: usize = 300_000_000; //real
// pub const MAX_VOCAB: usize = 20000000; // Debugging...

pub struct Dict {
    pub wordidx: Vec<AtomicCell<Option<usize>>>, // Randomish, based on hash
    // words: dashmap::DashMap<usize, Vec<u8>>,
    pub words: Vec<OnceCell<ThinVec<u8>>>,
    // words:   RwLock<Vec<Entry>>, // Sequential
    // words: Vec<AtomicCell<Option<Vec<u8>>>>,
    pub counts:  Vec<AtomicCell<usize>>, // Storing the counts separately now...
    pub tokens:  AtomicCell<usize>,
    pub size:    AtomicCell<u32>,
    pub entries: AtomicCell<usize>,
    // minimum_threshold: AtomicCell<u32>,
    // maximum_threshold: AtomicCell<u32>,
    // ignore_table: RwLock<Vec<bool>>,
    // pct75: u32,
    // discard_table: RwLock<Vec<usize>>
}

/*#[derive(Debug)]
pub struct Entry {
    pub kmer: Vec<u8>,
    pub count: AtomicCell<u32>,
}*/

impl Dict {
    pub fn new() -> Dict {
        let mut wordidx = Vec::new();
        wordidx.resize_with(MAX_VOCAB, || AtomicCell::new(None));

//        let rcidx = RwLock::new(Vec::new());
//        rcidx.write().unwrap().resize_with(MAX_VOCAB, || None);

        let mut counts = Vec::new();
        counts.resize_with(MAX_VOCAB, || AtomicCell::new(0));

        // let mut words = Vec::new();
        // words.resize_with(MAX_VOCAB, || AtomicCell::new(None));

        let mut words = Vec::with_capacity(MAX_VOCAB);
        words.resize_with(MAX_VOCAB, OnceCell::new);

        Dict {
            // vec![None; MAX_VOCAB],
            wordidx,
            // words: RwLock::new(Vec::with_capacity((MAX_VOCAB as f64 * 0.75) as usize)),
            // words: dashmap::DashMap::with_capacity(18, MAX_VOCAB),
            words,
            counts,
            // ignore_table: RwLock::new(vec![false; MAX_VOCAB]),
            tokens: AtomicCell::new(0),
            size: AtomicCell::new(0),
            entries: AtomicCell::new(0),
            // minimum_threshold: AtomicCell::new(1),
            // maximum_threshold: AtomicCell::new(std::u16::MAX as u32),
            // pct75: (0.75 * MAX_VOCAB as f64) as u32,
            // discard_table: RwLock::new(Vec::new())
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
        let mut rc: ThinVec<u8> = ThinVec::new();
        rc.extend_from_slice(&kmer);
        super::complement_nucleotides(&mut rc);
        rc.reverse();
        rc

    }

    #[inline]
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
        let x = wyhash(&kmer, 43_988_123) as usize % MAX_VOCAB;

        // let mut hasher = WyHash::with_seed(3);
        // hasher.write(&rc);
        // let y = hasher.finish() as usize % MAX_VOCAB;
        let y = wyhash(&rc, 43_988_123) as usize % MAX_VOCAB;

        // Very slow...
        // let x: usize = self.transmute_hash(&kmer) % MAX_VOCAB;
        // let y: usize = self.transmute_hash(&rc) % MAX_VOCAB;

        std::cmp::min(x,y)
    }

/*     fn transmute_hash(&self, data: &[u8]) -> usize {
        let mut x_raw = data[..].chunks_exact(8);
        let mut xrem = x_raw.remainder().to_vec();
        xrem.resize(8, 0);
        let mut x: u64 = unsafe { mem::transmute(&xrem) };
        x_raw.for_each(|xs| {
            unsafe { 
                x = x.wrapping_add(mem::transmute_copy::<_, u64>(&xs));
            };
        });

        x as usize
    } */

    /* fn get_id_prev(&self, kmer: &[u8]) -> usize {
        let mut id = (seahash::hash(kmer) as usize) % MAX_VOCAB;
        while self.wordidx[id] != None && self.words[self.wordidx[id].unwrap()].kmer != kmer {
            id = (id + 1) % MAX_VOCAB;
        }

        id as usize
    } */

    //pub fn add(&mut self, kmer: &[u8]) -> usize {
    pub fn add(&self, kmer: &[u8]) {
        self.tokens.fetch_add(1);

//        if *self.ignore_table.read().unwrap().get(hash).unwrap() {
//            return ()
//        }

        let id = self.get_id(kmer);

        if self.wordidx[id].load() == None {
            // self.entrieself.wordidxs += 2;

            // For DashMap
            // self.words.insert(id, kmer.to_vec());
            let mut word = ThinVec::with_capacity(kmer.len());
            word.extend_from_slice(&kmer);
            if let Err(val) = self.words[id].set(word) {
                self.add(&val); return; 
            }

            // match self.words.try_get_mut()

            // let mut words_w = self.words.write().unwrap();
            
            /* let wordidx = match self.discard_table.write().unwrap().pop() {
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
            drop(words_w); */

            /* match self.words[id].compare_and_swap(None, kmer) {
                None => (),
                Some => // Collision
                 { self.add(kmer);
                   return();
                 }
            }; */

            let wordidx = Some(self.size.fetch_add(1) as usize);
            
            match self.wordidx[id].compare_and_swap(None, wordidx) {
                None  => (),
                Some(_) => // Collision
                {
                    // So nevermind, start over and add this entry to the discard table...
                    // self.discard_table.write().unwrap().push(x);
                    self.add(kmer);
                    return;
                }
            }; // Points to word

            self.entries.fetch_add(1);


            // RC's are either +1 or -1 so we can't simply add one...

            // k of 13 gives LOTS of singletons...
            /* if self.size.load() >= self.pct75 {
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
                    self.ignore_table.write().unwrap()[self.calc_hash(&v.kmer, &self.get_rc(&v.kmer))] = true;
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
                    self.ignore_table.write().unwrap()[self.calc_hash(&v.kmer, &self.get_rc(&v.kmer))] = true;
                    self.discard_table.write().unwrap().push(k);
                }

                self.entries.fetch_sub(removed);

                println!("");
                println!("Done with the pruning...");
                println!("{}", self.size.load());
                println!("");
            } */

        }
        
        self.counts[self.wordidx[id].load().unwrap()].fetch_add(1);
/*
        let mut added = false;
        let backoff = Backoff::new();
        while !added {
            let words_r = self.words.try_lock();
            match words_r {
                Ok(words_r) => {
                    words_r[self.wordidx[id].load().unwrap()].count.fetch_add(1);
                    added = true;
                },
                _ => {}
            }

            if added {
                break;
            }

            backoff.spin(); */
           

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
        //}

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