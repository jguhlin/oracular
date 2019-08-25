use seahash;
// use std::sync::RwLock;

const MAX_VOCAB: usize = 500000000;

pub struct Dict {
    wordidx: Vec<Option<usize>>, // Randomish, based on hash
    rcidx:   Vec<Option<usize>>, // Given a standard word ID, get the ID of the reverse complement
    words:   Vec<Entry>, // Sequential
    pub tokens:  usize,
    pub size:    u32,

}

#[derive(Debug)]
pub struct Entry {
    pub kmer: Vec<u8>,
    pub count: u64,
}

impl Dict {
    pub fn new() -> Dict {
        Dict {
            wordidx: vec![None; MAX_VOCAB],
            rcidx: vec![None; MAX_VOCAB],
            words: Vec::with_capacity(MAX_VOCAB / 1000000),
            tokens: 0,
            size: 0,
        }
    }

    #[inline]
    fn get_id(&self, kmer: &[u8]) -> usize {
        let mut id = (seahash::hash(kmer) as usize) % MAX_VOCAB;
        while (self.wordidx[id] != None && self.words[self.wordidx[id].unwrap()].kmer != kmer) {
            id = (id + 1) % MAX_VOCAB;
        }

        id as usize
    }

    //pub fn add(&mut self, kmer: &[u8]) -> usize {
    pub fn add(&mut self, kmer: &[u8]) {
        self.tokens += 2;
        let id = self.get_id(&kmer);
        if (self.wordidx[id] == None) {
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

            if (self.wordidx[rcid] == None) {
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
    }

    /*
    pub fn add_to_id(&mut self, id: usize) {
        self.words[id].count += 1;
    } */

}