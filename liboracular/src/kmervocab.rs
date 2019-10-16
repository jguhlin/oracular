use serde::Serialize;

use std::collections::HashMap;

use finalfrontier::CountedType;
use finalfrontier::Vocab;

use std::hash::Hash;
use std::borrow::Borrow;

use finalfrontier::idx::{WordWithSubwordsIdx};
use finalfrontier::idx::WordIdx;

use finalfrontier::{NGramConfig};

use finalfusion::subword::{
    Indexer, NGramIndexer, SubwordIndices
};

use finalfusion::chunks::vocab::{VocabWrap, FinalfusionNGramVocab};

pub type Kmer = CountedType<String>;

/*
impl Kmer {
    pub fn new(label: String, count: usize) -> Self {
        Kmer { label, count }
    }
    pub fn count(&self) -> usize {
        self.count
    }
    pub fn label(&self) -> &String {
        &self.label
    }
} */


/// Hyperparameters for Subword vocabs.
#[derive(Clone, Copy, Debug, Serialize)]
#[serde(rename = "KmerVocab")]
#[serde(tag = "type")]
pub struct KmerVocabConfig<V> {
    pub min_count: u32,
    pub discard_threshold: f32,
    pub min_n: u32,
    pub max_n: u32,
    pub indexer: V,
}

#[derive(Clone)]
pub struct KmerVocab<C, I> {
    config: KmerVocabConfig<C>,
    words: Vec<Kmer>, // Switch to Vec<Kmer> ?
    indexer: I,
    subwords: Vec<Vec<u64>>,
    discards: Vec<f32>,
    index: HashMap<String, usize>,
    n_tokens: usize,
}

impl From<KmerVocab<NGramConfig, NGramIndexer>> for VocabWrap {
    fn from(v: KmerVocab<NGramConfig, NGramIndexer>) -> Self {
        let x = FinalfusionNGramVocab {
            indexer: v.indexer,
            indices: v.index,
            words: v.words.into_iter().map(|x| x.label().to_string()).collect::<Vec<String>>(),
            min_n: 9,
            max_n: 9
        };
        VocabWrap::FinalfusionNGramVocab(x)
    }
}

impl<C, I> KmerVocab<C, I>
where
    C: Copy + Clone,
    I: Indexer,
{
    pub fn new(
        config: KmerVocabConfig<C>,
        words: Vec<Kmer>,
        n_tokens: usize,
        indexer: I,
    ) -> Self {
        let index = create_indices(&words);
        let subwords = Self::create_subword_indices(
            config.min_n as usize,
            config.max_n as usize,
            &indexer,
            &words,
        );
        let discards = create_discards(config.discard_threshold, &words, n_tokens);
        KmerVocab {
            config,
            discards,
            indexer,
            words,
            subwords,
            index,
            n_tokens,
        }
    }

    fn create_subword_indices(
        min_n: usize,
        max_n: usize,
        indexer: &I,
        words: &[Kmer],
    ) -> Vec<Vec<u64>> {
        let mut subword_indices = Vec::new();

        for word in words {
            subword_indices.push(
                word.word()
                    .subword_indices(min_n, max_n, indexer)
                    .into_iter()
                    .map(|idx| idx + words.len() as u64)
                    .collect(),
            );
        }

        assert_eq!(words.len(), subword_indices.len());

        subword_indices
    }

    pub fn word(&self, word: &str) -> Option<&Kmer> {
        self.idx(word)
            .map(|idx| &self.words[idx.word_idx() as usize])
    }
}

impl<C, I> KmerVocab<C, I> {
    pub(crate) fn subword_indices_idx(&self, idx: usize) -> Option<&[u64]> {
        self.subwords.get(idx).map(|v| v.as_slice())
    }
}

impl<C, I> Vocab for KmerVocab<C, I>
where
    C: Copy + Clone,
    I: Indexer,
{
    type VocabType = String;
    type IdxType = WordWithSubwordsIdx;
    type Config = KmerVocabConfig<C>;

    fn config(&self) -> KmerVocabConfig<C> {
        self.config
    }

    fn idx<Q>(&self, key: &Q) -> Option<Self::IdxType>
    where
        Self::VocabType: Borrow<Q>,
        Q: Hash + ?Sized + Eq,
    {
        self.index.get(key).and_then(|idx| {
            self.subword_indices_idx(*idx)
                .map(|v| WordWithSubwordsIdx::new(*idx as u64, v))
        })
    }

    fn discard(&self, idx: usize) -> f32 {
        self.discards[idx]
    }

    fn n_input_types(&self) -> usize {
        self.len() + self.indexer.upper_bound() as usize
    }

    fn types(&self) -> &[Kmer] {
        &self.words
    }

    fn n_types(&self) -> usize {
        self.n_tokens
    }
}

fn create_indices(types: &[Kmer]) -> HashMap<String, usize>
{
    let mut token_indices = HashMap::new();

    for (idx, item) in types.iter().enumerate() {
        token_indices.insert(item.label().clone(), idx);
    }

    assert_eq!(types.len(), token_indices.len());

    token_indices
}

fn create_discards<S>(
    discard_threshold: f32,
    types: &[CountedType<S>],
    n_tokens: usize,
) -> Vec<f32> {
    let mut discards = Vec::with_capacity(types.len());

    for item in types {
        let p = item.count() as f32 / n_tokens as f32;
        let p_discard = discard_threshold / p + (discard_threshold / p).sqrt();

        // Not a proper probability, upper bound at 1.0.
        discards.push(1f32.min(p_discard));
    }

    discards
}