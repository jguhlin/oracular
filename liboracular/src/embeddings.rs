use finalfrontier::{
    Vocab, SkipgramTrainer, SGD, ModelType, LossType, SkipGramConfig, CommonConfig
};

use finalfusion::prelude::VocabWrap;

use rand_xorshift::XorShiftRng;
use rand::{SeedableRng, Rng};
use std::fs::File;
use std::io::{BufWriter};
use std::thread;
use std::sync::{Arc};
use crossbeam::atomic::AtomicCell;
use std::time::Duration;

use stdinout::OrExit;

use crossbeam::queue::{ArrayQueue};
use crossbeam::utils::Backoff;

use serde::Serialize;

use crate::threads::{sequence_generator, Sequence, ThreadCommand};
use opinionated::fasta::{complement_nucleotides};
use finalfrontier::WriteModelBinary;
// use finalfrontier::WriteModelText;
use finalfrontier::io::TrainInfo;

use indicatif::ProgressBar;
use indicatif::ProgressStyle;

use crossbeam::queue::{PushError};

pub fn train<V>(num_threads: usize, dims: usize, epochs: usize, context_size: usize, vocab: V, filename: &str, kmer_size: usize) -> ()
where
    V: Vocab<VocabType = String> + Into<VocabWrap> + Clone + Send + Sync + 'static,
    V::Config: Serialize,
    V: Vocab + Into<VocabWrap>,
    for<'a> &'a V::IdxType: IntoIterator<Item = u64>,
{
    println!("Creating embeddings...");

    let train_info = TrainInfo::new(
        filename.to_string(),
        "embeddings.embed".to_string(), // OUTPUT filename, make variable TODO
        num_threads as usize,);

    let skipgram_config = SkipGramConfig {
            context_size: context_size as u32, // Should be 6 (trying out 5 though)
            model: ModelType::SkipGram,            // loss: 0.12665372
            // model: ModelType::DirectionalSkipgram,// loss: 0.12277688
            // model: ModelType::StructuredSkipGram, // loss: 0.1216571

            // Loss above is lr 0.08, epochs: 10, dims: 32, neg_samples 100
            // For vvulg genome, context of 6

            // StructuredSkipGram -- Context Size of 10: 0.16051194
            //                       54:32.10elapsed 6014%CPU
        };

    let common_config = CommonConfig {
        loss: LossType::LogisticNegativeSampling,
        dims: dims as u32,
        epochs: epochs as u32,
        lr: 0.05, // Should be 0.07 for small datasets -- 0.12 for nt?
        negative_samples: 25, // Should be 20 - 100
        zipf_exponent: 0.05,
    };

    let trainer = SkipgramTrainer::new(
        vocab,
        XorShiftRng::from_entropy(),
        common_config,
        skipgram_config,
    );

    let sgd = SGD::new(trainer.into());
    let filename = filename.to_string();

    let msg = Arc::new("".to_string());

    let (seq_queue, jobs, generator_done, generator, mut children) 
        = sequence_generator(kmer_size, &filename, Arc::clone(&msg), common_config.epochs as u64);

    println!("Starting embedding worker threads...");

    for _ in 0..num_threads {
        let sgd = sgd.clone();
        let seq_queue = Arc::clone(&seq_queue);
        let jobs = Arc::clone(&jobs);
        let msg = Arc::clone(&msg);

        children.push(thread::spawn(move || {
            embedding_worker(
                sgd,
                common_config.lr,
                kmer_size,
                seq_queue,
                jobs,
                common_config.epochs as usize,
                msg,
            );
        }));
    }

    println!("Waiting for sequence generator to finish...");
    let backoff = Backoff::new();
    while !*generator_done.read().unwrap() {
        backoff.snooze();
    }

    generator.join().expect("Unable to join generator thread...");

    let update_interval = Duration::from_millis(500);

    backoff.snooze();

    let n_tokens = sgd.model().input_vocab().n_types();
    let pb = ProgressBar::new(n_tokens as u64 * common_config.epochs as u64);
    pb.set_style(ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:50.yellow} {eta_precise} {msg}")
                .progress_chars("█▇▆▅▄▃▂▁  "));
    
    // pb.enable_steady_tick(750);

    while !seq_queue.is_empty() {
        thread::sleep(update_interval);
        pb.set_position(sgd.n_tokens_processed() as u64);
        let lr = (1.0 - (sgd.n_tokens_processed() as f32 / (common_config.epochs as usize * n_tokens) as f32)) * common_config.lr;

        // pb.set_position(pb_len.saturating_sub(seq_queue.len() as u64));
        pb.set_message(&format!("loss/lr: {} {}", sgd.train_loss(), lr));
    }

    while jobs.load() > 0 {
        thread::sleep(update_interval);
        pb.set_position(sgd.n_tokens_processed() as u64);
        let lr = (1.0 - (sgd.n_tokens_processed() as f32 / (common_config.epochs as usize * n_tokens) as f32)) * common_config.lr;

        // pb.set_position(pb_len.saturating_sub(seq_queue.len() as u64));
        pb.set_message(&format!("loss/lr: {} {}", sgd.train_loss(), lr));
    }

    for _ in 0..num_threads {
	let mut i = 0;
        let mut result = seq_queue.push(ThreadCommand::Terminate);
        while let Err(PushError(wp)) = result {
		i = i + 1;
        	result = seq_queue.push(wp);
		thread::sleep(update_interval);
		if i > 1_000_000 {
            println!("Too many, giving up on issuing Terminate commands");
			break;
		}
        }

//        match seq_queue.push(ThreadCommand::Terminate) {
//            Ok(_) => (),
//            Err(x) => panic!("Unable to send command... {:#?}", x)
//        }
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

    println!();
    println!();
    println!("Final loss: {}", sgd.train_loss());

    // Wait until all threads have finished.

    // sgd.into_model();
        // .write_model_binary(&mut output_writer, train_info)
        // .or_exit("Cannot write model", 1);

    let mut output_writer = BufWriter::new(
        File::create("embeddings.embed").or_exit("Cannot open output file for writing.", 1),
    );

    println!("Tokens processed: {}", sgd.n_tokens_processed());

    thread::sleep(update_interval);

    let model = sgd.into_model();
    
    // let (x, y) = model.into_parts().unwrap();
    // Weird error now: thread 'main' panicked at 'called `Result::unwrap()` on an `Err` value: ErrorMessage { msg: "Cannot unwrap input matrix." }', src/libcore/result.rs:1165:

    // println!("{}", y);

    model.write_model_binary(&mut output_writer, train_info);

    // model

}

fn embedding_worker<R, V>(
    mut sgd: SGD<SkipgramTrainer<R, V>>,
    start_lr: f32,
    kmer_size: usize,
    seq_queue: Arc<ArrayQueue<ThreadCommand<Sequence>>>,
    jobs: Arc<AtomicCell<usize>>,
    epochs: usize,
    mut msg: Arc<String>
) where
    R: Clone + Rng,
    V: Vocab<VocabType = String>,
    V::Config: Serialize,
    for<'a> &'a V::IdxType: IntoIterator<Item = u64>,
{
    let backoff = Backoff::new();
    let mut kmer_steps: Vec<usize> = Vec::new();

    kmer_steps.push(0);
    let mut rng = rand::thread_rng();
    kmer_steps.push(rng.gen_range(1, kmer_size));
    let half_kmer = (kmer_size as f64 / 2.0).floor() as usize;

    kmer_steps.push(half_kmer);
    kmer_steps.push(kmer_size - 1);
    
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

            // TODO: Change lr curve function...
            let lr = (1.0 - (sgd.n_tokens_processed() as f32 / (epochs as usize * n_tokens) as f32)) * start_lr;

            let mut sentence: Vec<String> = Vec::new();

            rawseq.chunks_exact(kmer_size).for_each(|x| { 
                unsafe { sentence.push(String::from_utf8_unchecked(x.to_vec())); }
            });
            sgd.update_sentence(&sentence, lr);
            sentence.clear();

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

            match Arc::get_mut(&mut msg) {
                Some(x) => *x = format!("loss: {} lr: {}", sgd.train_loss(), lr),
                None    => ()
            };

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
