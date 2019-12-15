use finalfusion::prelude::*;
// use finalfusion::similarity::EmbeddingSimilarity;
use ndarray;
use snap;
use bincode;
use std::str;

use crate::kmer_counting::FinalDict;

use std::fs::File;
use std::io::BufReader;

use ndarray::{Array};

pub type Embedding = ndarray::ArrayBase<ndarray::OwnedRepr<f32>, ndarray::Dim<[usize; 1]>>;

pub struct EmbeddingResult {
    pub embedding: Embedding,
    pub count: usize
}

pub struct Model {
    vocab: FinalDict,
    pub embeddings: Embeddings<VocabWrap, StorageViewWrap>,
    pub k: usize
}

impl Model {

    pub fn new(model_filename: String, vocab_filename: String, k: usize) -> Model {
        let kmercounts_fh = snap::Reader::new(BufReader::with_capacity(1024 * 1024 * 1, File::open(vocab_filename.clone()).unwrap()));
        let mut reader = BufReader::new(File::open(model_filename).unwrap());

        Model {
            vocab: bincode::deserialize_from(&mut BufReader::new(kmercounts_fh)).expect("Unable to read to bincode file"),
            embeddings: Embeddings::mmap_embeddings(&mut reader).unwrap(),
            k: k as usize
        }
    }

    #[inline]
    pub fn get_embedding_and_count(&self, kmer: &str) -> EmbeddingResult {
        return EmbeddingResult { embedding: self.get_embedding(kmer), count: self.get_count(kmer) }
    }

    #[inline]
    pub fn get_count(&self, kmer: &str) -> usize {
        match self.vocab.words.get(kmer) {
            Some(x) => x.clone() as usize,
            None    => 0
        }
    }

    // Default to 1 instead of 0, for kmers not found in our database but used to create a weighting...
    // Set to 2 now to avoid 1/log10(1) being infinite...
    #[inline]
    pub fn get_count1(&self, kmer: &str) -> usize {
        match self.vocab.words.get(kmer) {
            Some(x) => x.clone() as usize + 2_usize,
            None    => 2
        }
    }

/*    #[inline]
    pub fn get_embedding_previous(&self, kmer: &str) -> Embedding {
           match self.embeddings.embedding(&kmer) {
            Some(x) => x.into_owned(),
            None    => panic!("Unable to get embeddings for {}", kmer)
        }
    } */

    #[inline]
    pub fn get_embedding(&self, kmer: &str) -> Embedding {
           match self.embeddings.embedding_with_norm(&kmer) {
            Some(x) => x.embedding.into_owned(),
            None    => panic!("Unable to get embeddings for {}", kmer)
            // Probably have to use the ngrams module to get indices, and then average...
        }
    }

    pub fn get_weighted_embedding(&self, seq: &str) -> Embedding {

        let seq = seq.to_ascii_uppercase();

        let mut counts: Vec<f64> = Vec::with_capacity(seq.len());
        let mut embeddings: Vec<Embedding> = Vec::with_capacity(seq.len());

        // This probably will slow things down, but does increase accuracy a little
        // 0.99307024 to 0.9999744
        // For a sequence where a single base has been removed...
        // for i in 0..self.vocab.k {
            // let i = i as usize;
            let mut end_of_seq_kmer = seq[seq.len() - self.vocab.k as usize..].to_string();

            // let chunks = seq[i..].as_bytes().chunks_exact(self.vocab.k as usize);
            let chunks = seq[..].as_bytes().chunks_exact(self.vocab.k as usize);
            let remainder = chunks.remainder();

            // Our ngrams are 9, so we can generate an embedding if it's at least 9 characters
            if remainder.len() >= 11 {
                end_of_seq_kmer = unsafe { String::from_utf8_unchecked(remainder.to_vec()) };
            } // Otherwise we use the last 11 and let it even out...

            for kmer in chunks {
                let kmer = unsafe { str::from_utf8_unchecked(kmer) };
                embeddings.push(self.get_embedding(kmer));
                counts.push(self.get_count1(kmer) as f64);
            }

            if remainder.len() > 0 {
                embeddings.push(self.get_embedding(&end_of_seq_kmer));
                counts.push(self.get_count1(&end_of_seq_kmer) as f64);
            }
        // }

        println!("Embeddings: {:#?}", embeddings);

        println!("Counts: {:#?}", counts);

        // NO: Values come out too evenly, magnitude is too large...
        // counts = counts.iter().map(|x| 1_f64 - (x/self.vocab.tokens as f64)).collect();

        let max = counts.iter().fold(0.0, |acc, &y| f64::max(acc, y));
        let sum: f64 = counts.iter().sum();
        // Not convinced this is the best, but it's "working" for now and should be tested...
        /* After softmax:
        [
            0.10703707332218637,
            0.10650027009741918,
            0.10602246164531375,
            0.10684551937196488,
            0.10693002502432691,
            0.03938646845697209,
            0.10658055279671563,
            0.10690757169372662,
            0.106918137379581,
            0.1068719202117935,
        ] */
        
        // counts = counts.iter().map(|x| 1_f64 - (x/max)).collect();
        
        // println!("{:#?}", counts);
        // Trying inverse of log10(freq)
        counts = counts.iter().map(|x| 1_f64/x.log10()).collect();
        
        /* Next two after softmax, so too much weight to a single kmer...
        [
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ] */
        // No to the following two.
        // counts = counts.iter().map(|x| (x - max).abs()).collect();
        // counts = counts.iter().map(|x| 1_f64/(x/max)).collect();

        println!("Counts: {:#?}", counts);
        let weights = softmax(counts);
        println!("Weights: {:#?}", weights);

        let mut embedding: Embedding = Array::zeros(embeddings[0].len());

        for (x,i) in embeddings.iter().zip(weights) {
            // println!("{:#?}", embedding);
            let l2_norm = euclidian_norm(&x);
            if l2_norm > 0.0 {
                embedding = embedding + ((x / l2_norm) * i as f32);
            }
        }

        /* for x in embeddings.iter() {
            let l2_norm = euclidian_norm(&x);
            if l2_norm > 0.0 {
                embedding = embedding + (x / l2_norm);
            }
        } */


        embedding
    }

    #[allow(dead_code)]
    pub fn cos_distance(a: &Embedding, b: &Embedding,) -> f32 {
        // let a = a.clone();
        // let b = b.clone();
        let upper = a.dot(b);
        let lower = euclidian_norm(a) * euclidian_norm(b);
        return upper / lower
    }

}

// Should have done using ndarray, but this is quick enough for now...
fn softmax(m: Vec<f64>) -> Vec<f64> {
    let max = m.iter().fold(0.0, |acc, &y| f64::max(acc, y));
    let y: Vec<f64> = m.iter().map(|x| (x - max).exp()).collect();
    let sum: f64 = y.iter().sum();
    y.iter().map(|x| x / sum).collect()
}

// Right now cos similarity is best, but I always seem to use this for something or other...
#[allow(dead_code)]
#[inline]
fn euclidian_distance(a: &Embedding, b: &Embedding,) -> f32 {
    let result = a - b;
    return result.dot(&result).sqrt();
}

#[allow(dead_code)]
fn euclidian_norm(a: &Embedding) -> f32 {
    return a.dot(a).sqrt();
}

// First, you missed the part that get_sentence_vector is not just a simple "average". Before FastText sum each word vector, each vector is divided with its norm (L2 norm) and then the averaging process only involves vectors that have positive L2 norm value.

// Second, a sentence always ends with an EOS. So if you try to calculate manually you need to put EOS before you calculate the average.

// try this (I assume the L2 norm of each word is positive):

/*
def l2_norm(x):
   return np.sqrt(np.sum(x**2))

def div_norm(x):
   norm_value = l2_norm(x)
   if norm_value > 0:
       return x * ( 1.0 / norm_value)
   else:
       return x

# Getting word vectors for 'one' and 'two'.
one = model.get_word_vector('yksi')
two = model.get_word_vector('kaksi')
eos = model.get_word_vector('\n')

# Getting the sentence vector for the sentence "one two" in Finnish.
one_two = model.get_sentence_vector('yksi kaksi')
one_two_avg = (div_norm(one) + div_norm(two) + div_norm(eos)) / 3
*/
