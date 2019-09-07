use crate::dictionary::Dict;
use finalfrontier::{SimpleVocab, SimpleVocabConfig, VocabBuilder};
use std::sync::{Arc};
use std::ops::Rem;

pub(crate) fn create_embeddings(dict: Arc<Dict>) {
    println!("Creating Vocab");
    let vocab: SimpleVocab<Vec<u8>> = build_vocab(dict);
}

fn build_vocab(dict: Arc<Dict>) -> SimpleVocab<Vec<u8>> {
    let mut builder: VocabBuilder<_, Vec<u8>> = VocabBuilder::new(
            SimpleVocabConfig { min_count: 1, discard_threshold: 0.00001 });

    let dict_iter = dict.words.iter().zip(dict.counts.iter());
    let mut total: usize = 0;
    for (word, count) in dict_iter {
        total += 1;
        if (total.checked_rem(1000000) == Some(0)) {
            println!("{}", total);
        }
        // TODO: Make it so we can pass an absolute value instead of looping...
        if let Some(wordu) = word.get() {
            for _ in 0..=count.load() {
                builder.count(wordu.to_vec());
            }
        }
    }
    
    builder.into()
    
}