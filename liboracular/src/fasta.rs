use crate::threads::{sequence_generator, Sequence, ThreadCommand, SequenceTargetContexts, SequenceBatch};

use crossbeam::queue::{ArrayQueue, PushError};
use std::sync::{Arc, RwLock};
use crossbeam::atomic::AtomicCell;
use opinionated::fasta::{complement_nucleotides};
use crossbeam::utils::Backoff;

use std::thread::JoinHandle;

use rand::thread_rng;
use rand::seq::SliceRandom;

use std::thread::Builder;

const WORKERSTACKSIZE: usize = 64 * 1024 * 1024;  // Stack size (needs to be > BUFSIZE + SEQBUFSIZE)

/// Meant to be called from pyracular, but can be called from other spots too
/// # parse_ntfasta_target_contexts
/// Parses a nucleotide FASTA file and its reverse complement into batches of 
/// the FASTA ID, a target kmer (in the middle of window_size kmers) and context kmers
/// meant to be used in the learning process.
/// 
/// So given a k=2 and window_size of 4
/// AC TG GC AC TG GC AC TG TT AA CC GG TT
/// 
/// Will result in TG as a target word, and AC, TG, GC, AC, GC, AC, TG, TT as context
/// words. Followed by the entire window increasing by 1 and repeating, as well as the reverse complement.
/// 
/// SequenceBatch type is Vec<SequenceTargetContexts>
/// 
/// So this function fills up a buffer with *buffer_size* sets of *batch_size* SequenceTargetContexts, waiting for
/// python or another function to pull from it.
/// 
/// THIS FOLLOWING IS A LIE
/// Shuffling is done on windows of *shuffle_buffer* size (Recommend: 200k or more)
/// The FASTA file should probably itself be pre-shuffled to help things out... (TEST THIS TODO)
/// 
pub fn parse_ntfasta_target_contexts(
        kmer_size: usize,
        filename: &str,
        window_size: usize,
        batch_size: usize,
        shuffle_buffer: usize,
        buffer_size: usize,
        num_threads: usize) 

        -> (Arc<ArrayQueue<ThreadCommand<SequenceBatch>>>,
            Arc<ArrayQueue<ThreadCommand<Sequence>>>,
            Arc<ArrayQueue<ThreadCommand<SequenceTargetContexts>>>,
            JoinHandle<()>,
            Arc<RwLock<bool>>,
            Arc<AtomicCell<usize>>,
            Vec<JoinHandle<()>>)
    
    {

    let unshuffled_queue   = Arc::new(ArrayQueue::<ThreadCommand<SequenceTargetContexts>>::new(shuffle_buffer));
    let batch_queue        = Arc::new(ArrayQueue::<ThreadCommand<SequenceBatch>>::new(buffer_size));

    let (seq_queue, jobs, generator_done, generator, mut children) 
        = sequence_generator(kmer_size, &filename, 1_u64);

    for _ in 0..num_threads {
        let seq_queue = Arc::clone(&seq_queue);
        let unshuffled_queue = Arc::clone(&unshuffled_queue);
        let kmer_size = kmer_size.clone();
        let window_size = window_size.clone();
        let jobs = Arc::clone(&jobs);

        let child = match Builder::new()
            .name("ntfasta_worker".into())
            .stack_size(WORKERSTACKSIZE)
            .spawn(move || ntfasta_worker_thread(kmer_size, window_size, seq_queue, unshuffled_queue, jobs)) {
                Ok(x)  => x,
                Err(y) => panic!("Unable to create ntfasta_worker_thread {}", y)
            };

        children.push(child);
    }

    // Thread to handle some minor shuffling and putting them into batches...

    {
        let batch_queue = Arc::clone(&batch_queue);
        let jobs = Arc::clone(&jobs);
        let unshuffled_queue = Arc::clone(&unshuffled_queue);

        let child = match Builder::new()
            .name("ShuffleAndBatch".into())
            .stack_size(WORKERSTACKSIZE)
            .spawn(move || shuffle_and_batch(batch_size, unshuffled_queue, batch_queue, jobs)) {
                Ok(x) => x,
                Err(y) => panic!("Unable to create shuffle and batch thread {}", y)
            };

        children.push(child);
    }

    // Need to return things (generator_done, generator, jobs, children) so that other threads/processes can handle them...

    (batch_queue, seq_queue, unshuffled_queue, generator, generator_done, jobs, children)
}

fn shuffle_and_batch(
    batch_size: usize, 
    unshuffled_queue: Arc<ArrayQueue<ThreadCommand<SequenceTargetContexts>>>,
    batch_queue: Arc<ArrayQueue<ThreadCommand<SequenceBatch>>>, 
    jobs: Arc<AtomicCell<usize>>) 
{

    let backoff = Backoff::new();

    let mut received: Vec<SequenceTargetContexts> = Vec::with_capacity(batch_size*10);

    loop {

        if let Ok(command) = unshuffled_queue.pop() {

            // We are finished, end the thread...
            if let ThreadCommand::Terminate = command {
                while received.len() > 0 {
                    received.shuffle(&mut thread_rng()); // Give it a "good" shuffle

                    let this_batch_size = 
                        if received.len() > batch_size {
                            batch_size
                        } else {
                            received.len()
                        };

                    let batch: Vec<SequenceTargetContexts>;
                    batch = received.drain(0..this_batch_size).collect();
    
                    let wp = ThreadCommand::Work(batch);
    
                    let mut result = batch_queue.push(wp);
                    while let Err(PushError(wp)) = result {
                        backoff.snooze();
                        result = batch_queue.push(wp);
                    }
                }
                
                return;
            }

            received.push(command.unwrap());

            jobs.fetch_sub(1); // Yes? No? Not keeping tracking of it right now...

            if received.len() >= batch_size*10 {
                received.shuffle(&mut thread_rng()); // Give it a "good" shuffle

                let batch: Vec<SequenceTargetContexts>;
                batch = received.drain(0..batch_size).collect();

                let wp = ThreadCommand::Work(batch);

                let mut result = batch_queue.push(wp);
                while let Err(PushError(wp)) = result {
                    backoff.snooze();
                    result = batch_queue.push(wp);
                }
            }

        } else {
            backoff.snooze();
            backoff.reset();
        }
    }
}


/// This function does the work of splitting up kmers into target and contexts...
/// and shuffling it off to raw_queue. Reverse complement is handled in io_worker!
fn ntfasta_worker_thread (
    kmer_size: usize, 
    window_size: usize,
    seq_queue: Arc<ArrayQueue<ThreadCommand<Sequence>>>,
    unshuffled_queue: Arc<ArrayQueue<ThreadCommand<SequenceTargetContexts>>>,
    jobs: Arc<AtomicCell<usize>>) {

    let backoff = Backoff::new();

    loop {
        if let Ok(command) = seq_queue.pop() {

            // We are finished, end the thread...
            if let ThreadCommand::Terminate = command {
                return;
            }          

            let packet = command.unwrap();
            let rawseq = packet.rawseq;
            let id     = packet.id;

            jobs.fetch_sub(1);

            let length = rawseq.len();

            if length <= kmer_size*window_size*2+kmer_size {
                continue
            }

            let query_length = ((window_size * 2) + 1) * kmer_size;
            let end          = length - query_length;

            for i in 0..end {
                let seq = rawseq[i..i+query_length].to_vec();

                let mut contexts: Vec<String> = Vec::with_capacity(window_size*2);
                let target  : String;

                let kmers: Vec<&[u8]> = seq.chunks_exact(kmer_size).collect();

                for z in 0..window_size {
                    contexts.push(String::from_utf8(kmers[z].to_vec()).expect("Invalid UTF-8 Encoding"));
                }

                for z in 0..window_size {
                    contexts.push(String::from_utf8(kmers[window_size+z].to_vec()).expect("Invalid UTF-8 Encoding"));
                }

                target = String::from_utf8(kmers[window_size].to_vec()).expect("Invalid UTF-8 Encoding");

                let seq = SequenceTargetContexts { 
                        target,
                        contexts,
                        id: id.clone(),};

                jobs.fetch_add(1 as usize);
                let wp = ThreadCommand::Work(seq);
    
                let mut result = unshuffled_queue.push(wp);
                while let Err(PushError(wp)) = result {
                    backoff.snooze();
                    result = unshuffled_queue.push(wp);
                }
            }

        } else {
            backoff.snooze();
            backoff.reset();
        }
    }
}
