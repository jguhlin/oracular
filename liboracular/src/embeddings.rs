use finalfrontier::{
    SentenceIterator, SimpleVocab, SkipgramTrainer, SubwordVocab, Vocab, VocabBuilder,
    WriteModelBinary, SGD, SubwordVocabConfig, NGramConfig, ModelType, LossType, SkipGramConfig,
    CommonConfig
};

use finalfusion::prelude::VocabWrap;

use rand_xorshift::XorShiftRng;
use rand::{FromEntropy, Rng};
use crate::rand::SeedableRng;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::thread;

use std::time::Duration;

use stdinout::OrExit;

use serde::Serialize;

const PROGRESS_UPDATE_INTERVAL: u64 = 200;

pub fn train<V>(vocab: V, filename: &str)
where
    V: Vocab<VocabType = String> + Into<VocabWrap> + Clone + Send + Sync + 'static,
    V::Config: Serialize,
    for<'a> &'a V::IdxType: IntoIterator<Item = u64>,
{
    let num_threads = 48;

    let jobs = Arc::new(AtomicCell::new(0 as usize));
    let seq_queue = Arc::new(ArrayQueue::<ThreadCommand<Sequence>>::new(1024 * 128));
    let rawseq_queue = Arc::new(ArrayQueue::<ThreadCommand<Sequence>>::new(2048));
    let done = Arc::new(RwLock::new(false));
    let generator_done = Arc::new(RwLock::new(false));

    let mut output_writer = BufWriter::new(
        File::create("embeddings").or_exit("Cannot open output file for writing.", 1),
    );

    let skipgram_config = SkipGramConfig {            
            context_size: 8,
            model: ModelType::SkipGram,
        };

    let common_config = CommonConfig {
        loss: LossType::LogisticNegativeSampling,
        dims: 32,
        epochs: 5,
        lr: 0.05,
        negative_samples: 15,
        zipf_exponent: 0.5,
    };

    let trainer = SkipgramTrainer::new(
        vocab,
        XorShiftRng::from_entropy(),
        common_config,
        skipgram_config,
    );

    let sgd = SGD::new(trainer.into());
    let filename = filename.to_string();

    let mut children = Vec::with_capacity(num_threads);
    for thread in 0..num_threads {
        let sgd = sgd.clone();

        children.push(thread::spawn(move || {
            do_work(
                sgd,
                thread,
                n_threads,
                common_config.lr,
            );
        }));
    }

    show_progress(
        &common_config,
        &sgd,
        Duration::from_millis(PROGRESS_UPDATE_INTERVAL),
    );

    // Wait until all threads have finished.
    for child in children {
        let _ = child.join();
    }

//    sgd.into_model()
//        .write_model_binary(&mut output_writer)
//        .or_exit("Cannot write model", 1);
}

fn do_work<R, V>(
    mut sgd: SGD<SkipgramTrainer<R, V>>,
    thread: usize,
    n_threads: usize,
    start_lr: f32,
) where
    R: Clone + Rng,
    V: Vocab<VocabType = String>,
    V::Config: Serialize,
    for<'a> &'a V::IdxType: IntoIterator<Item = u64>,
{
    let n_tokens = sgd.model().input_vocab().n_types();

    let f = File::open(corpus_path.into()).or_exit("Cannot open corpus for reading", 1);
    let (data, start) =
        thread_data_text(&f, thread, n_threads).or_exit("Could not get thread-specific data", 1);

    let mut sentences = SentenceIterator::new(&data[start..]);
    while sgd.n_tokens_processed() < epochs as usize * n_tokens {
        let sentence = if let Some(sentence) = sentences.next() {
            sentence
        } else {
            sentences = SentenceIterator::new(&*data);
            sentences
                .next()
                .or_exit("Iterator does not provide sentences", 1)
        }
        .or_exit("Cannot read sentence", 1);

        let lr = (1.0 - (sgd.n_tokens_processed() as f32 / (epochs as usize * n_tokens) as f32))
            * start_lr;

        sgd.update_sentence(&sentence, lr);
    }
}