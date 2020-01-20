use crate::threads::{sequence_generator, Sequence, ThreadCommand};

use crossbeam::queue::{ArrayQueue, PushError};
use std::sync::{Arc, RwLock};
use crossbeam::atomic::AtomicCell;
use opinionated::fasta::{complement_nucleotides};
use crossbeam::utils::Backoff;

use std::thread::Builder;

pub struct SequenceTargetContexts {
    pub id: String,
    pub target: Vec<u8>,
    pub contexts: Vec<Vec<u8>>
}

pub type SequenceBatch = Vec<SequenceTargetContexts>;

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
/// Shuffling is done on windows of 200k entries...
/// The FASTA file should probably itself be pre-shuffled to help things out... (TEST THIS TODO)
/// 
pub fn parse_ntfasta_target_contexts(
        kmer_size: usize,
        filename: &str,
        window_size: usize,
        batch_size: usize,
        buffer_size: usize,
        num_threads: usize) 

        -> Arc<ArrayQueue<ThreadCommand<SequenceBatch>>> {
    

    let raw_queue   = Arc::new(ArrayQueue::<ThreadCommand<SequenceTargetContexts>>::new(buffer_size));
    let batch_queue = Arc::new(ArrayQueue::<ThreadCommand<SequenceBatch>>::new(buffer_size));

    let (seq_queue, jobs, generator_done, generator, mut children) 
        = sequence_generator(kmer_size, &filename, Arc::new("".to_string()), 1_u64);

    for _ in 0..num_threads {
        let seq_queue = Arc::clone(&seq_queue);
        let raw_queue = Arc::clone(&raw_queue);
        let kmer_size = kmer_size.clone();
        let window_size = window_size.clone();
        let jobs = Arc::clone(&jobs);

        let child = match Builder::new()
            .name("ntfasta_worker".into())
            .stack_size(WORKERSTACKSIZE)
            .spawn(move || ntfasta_worker_thread(kmer_size, window_size, seq_queue, raw_queue, jobs)) {
                Ok(x)  => x,
                Err(y) => panic!("Unable to create ntfasta_worker_thread {}", y)
            };

        children.push(child);
    }

    batch_queue
}

/// This function does the work of splitting up kmers into target and contexts...
/// and shuffling it off to raw_queue
fn ntfasta_worker_thread (
    kmer_size: usize, 
    window_size: usize,
    seq_queue: Arc<ArrayQueue<ThreadCommand<Sequence>>>,
    raw_queue: Arc<ArrayQueue<ThreadCommand<SequenceTargetContexts>>>,
    jobs: Arc<AtomicCell<usize>>) {

    let backoff = Backoff::new();

    loop {
        if let Ok(command) = seq_queue.pop() {

            // We are finished, end the thread...
            if let ThreadCommand::Terminate = command {
                return;
            }          

            let rawseq = command.unwrap();
            jobs.fetch_sub(1);

            let length = rawseq.len();
            let end = length - ((window_size * 2 - 1) * kmer_size);

            for i in 0..kmer_size {
                let kmer_bunches = rawseq[i..].chunks_exact(kmer_size).take(window_size*2+1).collect::<Vec<&[u8]>>();
                for bunch in kmer_bunches {
                    let contexts = Vec::with_capacity(window_size*2);
                    contexts.extend(bunch[0..window_size]);
                    contexts.extend(bunch[window_size+1..]);

                    let seq = SequenceTargetContexts { 
                        target: bunch[window_size],
                        contexts,
                        id,
                     }
                }
            }

            // TODO: Reverse complement

        } else {
            backoff.snooze();
            backoff.reset();
        }
    }
}