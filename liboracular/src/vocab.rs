use finalfrontier::{
    NGramConfig
};

use crate::kmer_counting::FinalDict;
use finalfrontier::Word;
use finalfusion::subword::{
    ExplicitIndexer, NGrams,
};

use std::collections::HashMap;
use crate::kmervocab::{KmerVocabConfig, KmerVocab};

pub fn build_vocab_from_finaldict(dict: FinalDict, mincount: usize, min_n: usize, max_n: usize) -> KmerVocab<NGramConfig, ExplicitIndexer>

// where V: Vocab<VocabType = String> + From<VocabBuilder<SubwordVocabConfig<NGramConfig>, String>>,
{

    let config = KmerVocabConfig {
                discard_threshold: 1e-4,
                min_count: mincount as u32, // Min count 2 or 3 for small datasets, 10 for nt
                max_n: max_n as u32,
                min_n: min_n as u32,
                indexer: NGramConfig { min_ngram_count: 5 }, // 5 for small datasets, 50 for nt
    };

    let mut words: Vec<_> = Vec::with_capacity(dict.entries as usize);
    let mut ngram_counts: HashMap<String, usize> = HashMap::new();

    let mut i: usize = 0;

    for (word, count) in dict.words {
        if count >= config.min_count as u64 {
            i = i + 1;
/*            let word_; // As a string instead of [u8]
            unsafe { 
                word_ = std::str::from_utf8_unchecked(&word).to_string();
            } */
            for ngram in NGrams::new(&word, config.min_n as usize, config.max_n as usize)
                .map(|ngram| ngram.to_string())
            {
                let cnt = ngram_counts.entry(ngram).or_default();
                *cnt += count as usize;
            }
            words.push(Word::new(word, count as usize));
        }
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

    println!("Final vocab size: {}", i);

    println!("Vocab Use Ratio: {}", (i as f32 / dict.entries as f32));

    KmerVocab::new(
            config,
            words,
            dict.tokens as usize,
            ExplicitIndexer::new(ngrams))

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
