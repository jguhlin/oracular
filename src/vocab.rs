use finalfrontier::Vocab;
use thincollections::thin_vec::ThinVec;
use finalfrontier::CountedType;
use std::hash::Hash;
use crate::dictionary::Dict;
use std::iter::FromIterator;

pub struct SimpleKmerVocabConfig {
    pub min_count: u64,
    pub discard_threshold: f32,
}

#[derive(Hash, PartialEq, Eq)]
pub struct Kmer {
    kmer: ThinVec<u8>,
}

pub struct SimpleKmerVocab {
    config: SimpleKmerVocabConfig,
    types: ThinVec<CountedType<usize>>,
    dict: Dict,
    // wordidx: ThinVec<usize>,
    // words: ThinVec<ThinVec<u8>>,
    // counts: ThinVec<usize>,
    discards: ThinVec<f32>,
    // n_tokens: usize,
    // n_entries: usize,
}



// Nah, just use the underlying dict..
/* impl SimpleKmerVocab {
    pub fn new(config: SimpleKmerVocabConfig, dict: Dict) -> SimpleKmerVocab {

        let wordidx = ThinVec::from_iter(dict.wordidx.iter().map(|x| x.load().unwrap()));
        let words = ThinVec::from_iter(dict.words.iter().map(|x| x.get().unwrap().clone() ));
        let counts = ThinVec::from_iter(dict.counts.iter().map(|x| x.load()));

        SimpleKmerVocab {
            config,
            wordidx, 
            words, 
            counts, 
            n_types: 0, 
            discards: ThinVec::new(), 
            n_tokens: dict.tokens.load(), 
            n_entries: dict.entries.load(),
            
        }
    }

    // get a specific context?? do they mean word?
    // pub fn get

} */

impl Vocab for SimpleKmerVocab {
    type VocabType = Kmer;
    // type IdxType 
    type Config = SimpleKmerVocabConfig;
    
    fn config(&self) -> SimpleKmerVocabConfig {
        self.config
    }

    fn idx(&self, kmer: &[u8]) -> Option<u64> {
        self.dict.get_id(kmer) as u64
    }
    
}