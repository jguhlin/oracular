use crate::dictionary::Dict;
use finalfrontier::{CommonConfig, ModelType, LossType, SkipgramTrainer, SkipGramConfig, SimpleVocab, SimpleVocabConfig, VocabBuilder, Vocab};
use std::sync::{Arc};
use finalfusion::vocab::VocabWrap;
use std::ops::Rem;
use serde::Serialize;
use rand::{Rng, FromEntropy};
use rand::rngs::StdRng;
use std::str;


pub(crate) fn create_embeddings(dict: Arc<Dict>) {
    println!("Creating Vocab");
    let vocab: SimpleVocab<String> = build_vocab(dict);
    train(vocab);
}

fn build_vocab(dict: Arc<Dict>) -> SimpleVocab<String>  {
    let mut builder: VocabBuilder<_, String> = VocabBuilder::new(
            SimpleVocabConfig { min_count: 2, discard_threshold: 0.00001 });

    let dict_iter = dict.words.iter().zip(dict.counts.iter());
    let mut total: usize = 0;
    for (word, count) in dict_iter {
        total += 1;
        if total.checked_rem(10000000) == Some(0) {
            println!("{}", total);
        }
        // TODO: Make it so we can pass an absolute value instead of looping...
        if let Some(wordu) = word.get() {
            for _ in 0..=count.load() {
                builder.count(unsafe { str::from_utf8_unchecked(&wordu) });
            }
        }
    }
    
    builder.into()
}

fn train<V>(vocab: V)
where
    V: Vocab<VocabType = String> + Into<VocabWrap> + Clone + Send + Sync + 'static,
    // V::Config: Serialize,
    // for<'a> &'a V::IdxType: IntoIterator<Item = u64>,
{
    let n_threads = 32;

    let common_config = CommonConfig {
        loss: LossType::LogisticNegativeSampling,
        dims: 16,
        epochs: 4,
        lr: 0.05,
        negative_samples: 5,
        zipf_exponent: 0.5,
    };

    let skipgram_config = SkipGramConfig {
        context_size: 10,
        model: ModelType::StructuredSkipGram,
    };

    let trainer = SkipgramTrainer::new(vocab, 
                            StdRng::from_entropy(),
                            common_config,
                            skipgram_config);






}