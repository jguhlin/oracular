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

        -> Arc<ArrayQueue<ThreadCommand<SequenceBatch>>> {
    

    let unshuffled_queue   = Arc::new(ArrayQueue::<ThreadCommand<SequenceTargetContexts>>::new(buffer_size));
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

    // Thread to handle shuffling and submitting batches...

    // Need to return things (generator_done, generator, jobs, children) so that other threads/processes can handle them...

    batch_queue
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
            let query_length = ((window_size * 2) + 1) * kmer_size;
            let end          = length - query_length;

            for i in 0..end {
                let seq = rawseq[i..i+query_length].to_vec();

                let mut contexts: Vec<Vec<u8>> = Vec::with_capacity(window_size*2);
                let target  : Vec<u8>;

                let kmers: Vec<&[u8]> = seq.chunks_exact(kmer_size).collect();

                for z in 0..window_size {
                    contexts.push(kmers[z].to_vec());
                }

                for z in 0..window_size {
                    contexts.push(kmers[window_size+z].to_vec());
                }

                target = kmers[window_size].to_vec();

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