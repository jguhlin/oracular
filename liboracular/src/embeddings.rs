use finalfrontier::{
    Vocab, SkipgramTrainer, SGD, ModelType, LossType, SkipGramConfig, CommonConfig
};

use finalfusion::prelude::VocabWrap;

use rand_xorshift::XorShiftRng;
use rand::{FromEntropy, Rng};
use std::fs::File;
use std::io::{BufWriter};
use std::thread;
use std::sync::{Arc};
use crossbeam::atomic::AtomicCell;

use stdinout::OrExit;

use crossbeam::queue::{ArrayQueue, PushError};
use crossbeam::utils::Backoff;

use serde::Serialize;

use crate::threads::{sequence_generator, Sequence, ThreadCommand};
use opinionated::fasta::{complement_nucleotides, capitalize_nucleotides};

pub fn train<V>(vocab: V, filename: &str, kmer_size: usize)
where
    V: Vocab<VocabType = String> + Into<VocabWrap> + Clone + Send + Sync + 'static,
    V::Config: Serialize,
    for<'a> &'a V::IdxType: IntoIterator<Item = u64>,
{
    let num_threads = 48;

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

    let (seq_queue, jobs, generator_done, generator, mut children) 
        = sequence_generator(kmer_size, &filename);

    let mut children = Vec::with_capacity(num_threads);
    for _ in 0..num_threads {
        let sgd = sgd.clone();
        let seq_queue = Arc::clone(&seq_queue);
        let jobs = Arc::clone(&jobs);

        children.push(thread::spawn(move || {
            embedding_worker(
                sgd,
                common_config.lr,
                kmer_size,
                seq_queue,
                jobs,
                common_config.epochs as usize,
            );
        }));
    }

    let backoff = Backoff::new();
    while !*generator_done.read().unwrap() {
        backoff.snooze();
    }

    generator.join().expect("Unable to join generator thread...");

    while !seq_queue.is_empty() {
        backoff.snooze();
    }

    println!("Seq queue is empty, sending terminate command...");

    while jobs.load() > 0 {
        backoff.snooze();
        // println!("{}", jobs.load());
    }

    for _ in 0..num_threads {
        match seq_queue.push(ThreadCommand::Terminate) {
            Ok(_) => (),
            Err(x) => panic!("Unable to send command... {:#?}", x)
        }
    }

    println!("Terminate commands sent, joining worker threads");

    for child in children {
        match child.join() {
            Ok(_) => (),
            Err(x) => panic!("Error joining worker thread... {:#?}", x)
        }
    }

    println!("Worker threads joined, getting sgd stats...");

/*    show_progress(
        &common_config,
        &sgd,
        Duration::from_millis(PROGRESS_UPDATE_INTERVAL),
    ); */

    // Wait until all threads have finished.

//    sgd.into_model()
//        .write_model_binary(&mut output_writer)
//        .or_exit("Cannot write model", 1);
}

fn embedding_worker<R, V>(
    mut sgd: SGD<SkipgramTrainer<R, V>>,
    start_lr: f32,
    kmer_size: usize,
    seq_queue: Arc<ArrayQueue<ThreadCommand<Sequence>>>,
    jobs: Arc<AtomicCell<usize>>,
    epochs: usize,
) where
    R: Clone + Rng,
    V: Vocab<VocabType = String>,
    V::Config: Serialize,
    for<'a> &'a V::IdxType: IntoIterator<Item = u64>,
{
    let backoff = Backoff::new();
    let mut kmer_steps: Vec<usize> = Vec::new();

    kmer_steps.push(0);
    kmer_steps.push(1);
    let half_kmer = (kmer_size as f64 / 2.0).floor() as usize;
    kmer_steps.push(half_kmer);
    
    loop {
        if let Ok(command) = seq_queue.pop() {
            // We are finished, end the thread...
            if let ThreadCommand::Terminate = command {
                return;
            }

            let rawseq = command.unwrap();
            jobs.fetch_sub(1);

            // Generate the kmer "sentences" here...
            // Best to use an iterator, probably...

            let n_tokens = sgd.model().input_vocab().n_types();

            let lr = (1.0 - (sgd.n_tokens_processed() as f32 / (epochs as usize * n_tokens) as f32)) * start_lr;

            let mut sentence: Vec<String> = Vec::new();

            for epoch in 0..epochs {
                
                for i in kmer_steps.iter() {
                    rawseq[*i+epoch..].chunks_exact(kmer_size).for_each(|x| { 
                        unsafe { sentence.push(String::from_utf8_unchecked(x.to_vec())); }
                    });
                    sgd.update_sentence(&sentence, lr);
                    sentence.clear();
                }

                let mut rc = rawseq.clone();
                complement_nucleotides(&mut rc);
                rc.reverse();
                
                for i in kmer_steps.iter() {
                    rawseq[*i..].chunks_exact(kmer_size).for_each(|x| { 
                        unsafe { sentence.push(String::from_utf8_unchecked(x.to_vec())); }
                    });
                    sgd.update_sentence(&sentence, lr);
                    sentence.clear();
                }
            }
        } else {
            backoff.snooze();
            backoff.reset();
        }
    }

/*    let n_tokens = sgd.model().input_vocab().n_types();

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
     */
}