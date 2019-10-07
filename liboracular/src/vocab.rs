use finalfrontier::{
    SentenceIterator, SimpleVocab, SkipgramTrainer, SubwordVocab, Vocab, VocabBuilder,
    WriteModelBinary, SGD, SubwordVocabConfig, NGramConfig
};

use crate::kmer_counting::FinalDict;

pub fn build_vocab_from_finaldict<V>(dict: FinalDict) -> V
where V: Vocab<VocabType = String> + From<VocabBuilder<SubwordVocabConfig<NGramConfig>, String>>,
{

    let config = SubwordVocabConfig {
                discard_threshold: 1e-4,
                min_count: 5,
                max_n: 11,
                min_n: 9,
                indexer: NGramConfig { min_ngram_count: 50 },
    };

    let mut builder = VocabBuilder::new(config);

    let mut i = 0;
    let mut total = 0;

    for (word, count) in dict.words {
        i = i + 1;
        
        if (i % 1_000_000) == 0 {
            println!("{} / {}", i, total);
        }

        for _ in 0..count {
            total = total + 1;
            builder.count(std::str::from_utf8(&word).expect("Kmer is an invalid string"));
        }

    }

    builder.into()

}
