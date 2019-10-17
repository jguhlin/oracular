use finalfrontier::{
    NGramConfig
};

use crate::kmer_counting::FinalDict;
use finalfrontier::Word;
use finalfusion::subword::{
    NGramIndexer, NGrams,
};

use std::collections::HashMap;
use crate::kmervocab::{KmerVocabConfig, KmerVocab};
use finalfrontier::WriteModelBinary;

pub fn build_vocab_from_finaldict(dict: FinalDict) -> KmerVocab<NGramConfig, NGramIndexer>

// where V: Vocab<VocabType = String> + From<VocabBuilder<SubwordVocabConfig<NGramConfig>, String>>,
{

    let config = KmerVocabConfig {
                discard_threshold: 1e-4,
                min_count: 5,
                max_n: 9,
                min_n: 9,
                indexer: NGramConfig { min_ngram_count: 5 },
    };

    let mut words: Vec<_> = Vec::with_capacity(dict.size as usize);
    let mut ngram_counts: HashMap<String, usize> = HashMap::new();

    for (word, count) in dict.words {
        let word_; // As a string instead of [u8]
        unsafe { 
            word_ = std::str::from_utf8_unchecked(&word).to_string();
        }
        for ngram in NGrams::new(&word_, config.min_n as usize, config.max_n as usize)
            .map(|ngram| ngram.to_string())
        {
            let cnt = ngram_counts.entry(ngram).or_default();
            *cnt += count as usize;
        }
        words.push(Word::new(word_, count as usize));
    }


    words.sort_unstable_by(|w1, w2| w2.cmp(&w1));

    let mut ngrams = ngram_counts
        .iter()
        .filter(|(_, count)| **count >= config.indexer.min_ngram_count as usize)
        .map(|(ngram, _)| ngram.clone())
        .collect::<Vec<_>>();

    ngrams.sort_unstable_by(|ngram1, ngram2| {
        let ngram1_cnt = ngram_counts[ngram1];
        let ngram2_cnt = ngram_counts[ngram2];
        (ngram2_cnt, ngram2).cmp(&(ngram1_cnt, ngram1))
    });

    KmerVocab::new(
            config,
            words,
            dict.size as usize,
            NGramIndexer::new(ngrams))

/*    let mut builder = VocabBuilder::new(config);

    let mut i = 0;
    let mut total: usize = 0;

       
    for (word, count) in dict.words {
        i = i + 1;
        
        if (i % 1_000_000) == 0 {
            println!("{} / {}", i, total);
        }

        for _ in 0..count {
            total = total.saturating_add(1);
            builder.count(std::str::from_utf8(&word).expect("Kmer is an invalid string"));
        }

    } 

    builder.into() */



}