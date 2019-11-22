use crossbeam::atomic::AtomicCell;
use opinionated::fasta::{capitalize_nucleotides};

use std::sync::{Arc, RwLock};

use std::thread;
use std::thread::Builder;
use std::thread::JoinHandle;

use std::fs::File;
use std::io::{BufReader, Read, BufRead};

use crossbeam::queue::{ArrayQueue, PushError};
use crossbeam::utils::Backoff;

use indicatif::ProgressBar;
use indicatif::ProgressStyle;

const STACKSIZE: usize = 256 * 1024 * 1024;  // Stack size (needs to be > BUFSIZE + SEQBUFSIZE)
const WORKERSTACKSIZE: usize = 64 * 1024 * 1024;  // Stack size (needs to be > BUFSIZE + SEQBUFSIZE)

pub type Sequence = Vec<u8>;

pub enum ThreadCommand<T> {
    Work(T),
    Terminate,
}

impl ThreadCommand<Sequence> {
    // Consumes the ThreadCommand, which is just fine...
    pub fn unwrap(self) -> Sequence {
        match self {
            ThreadCommand::Work(x)   => x,
            ThreadCommand::Terminate => panic!("Unable to unwrap terminate command"),
        }
    }
}

// Takes a file and submits Sequence type to buffers...
// Fills up the seq_buffer that it returns
// UP TO the calling fn/thread to handle clean-up...
pub fn sequence_generator(
    kmer_size: usize,
    filename: &str,
    msg: Arc<String>,
    epochs: u64
) -> (Arc<ArrayQueue<ThreadCommand<Sequence>>>, Arc<AtomicCell<usize>>, Arc<RwLock<bool>>, JoinHandle<()>, Vec<JoinHandle<()>>)
// seq_queue jobs generator_done generator children

{
    let jobs = Arc::new(AtomicCell::new(0 as usize));
    let seq_queue = Arc::new(ArrayQueue::<ThreadCommand<Sequence>>::new(256));
    let rawseq_queue = Arc::new(ArrayQueue::<ThreadCommand<Sequence>>::new(256));
    let generator_done = Arc::new(RwLock::new(false));

    let generator;

    let mut children = Vec::new();

    let filename = filename.to_string();

    // IO-bound, so 4 threads never seems to hurt anything (but seems to help)
    for _ in 0..4 {
        let seq_queue = Arc::clone(&seq_queue);
        let rawseq_queue = Arc::clone(&rawseq_queue);
        let jobs = Arc::clone(&jobs);

        let child = match Builder::new()
                        .name("IOWorker".into())
                        .stack_size(WORKERSTACKSIZE)
                        .spawn(move || io_worker_thread(kmer_size, rawseq_queue, seq_queue, jobs)) {
                            Ok(x)  => x,
                            Err(y) => panic!("{}", y)
                        };
        
        children.push(child);
    }

    { // Explicit lifetime
        let generator_done = Arc::clone(&generator_done);
        let rawseq_queue = Arc::clone(&rawseq_queue);
        let jobs = Arc::clone(&jobs);
        let mut buffer: Sequence = Vec::with_capacity(1024);

        generator = thread::Builder::new()
                            .name("Generator".to_string())
                            .stack_size(STACKSIZE)
                            .spawn(move||
        {
            for epoch in 0..epochs {
                let mut seqbuffer: Sequence = Vec::with_capacity(8 * 1024 * 1024); // 8 Mb to start, will likely increase...
                let mut seqlen: usize = 0;

                let file = match File::open(&filename) {
                    Err(why) => panic!("Couldn't open {}: {}", filename, why.to_string()),
                    Ok(file) => file,
                };

                let pb = ProgressBar::new(file.metadata().unwrap().len());

                let file = BufReader::with_capacity(64 * 1024 * 1024, file);

                pb.set_style(ProgressStyle::default_bar()
                    .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>5}/{len:5} {eta_precise} {msg}")
                    .progress_chars("█▇▆▅▄▃▂▁  "));

                let fasta: Box<dyn Read> = if filename.ends_with("gz") {
                    Box::new(flate2::read::GzDecoder::new(pb.wrap_read(file)))
                } else {
                    Box::new(pb.wrap_read(file))
                };

                let mut reader = BufReader::with_capacity(128 * 1024 * 1024, fasta);

                let backoff = Backoff::new();

                while let Ok(bytes_read) = reader.read_until(b'\n', &mut buffer) {
                // pb.set_message("hello");
                    if bytes_read == 0 {
                        // File is empty, we are done!
                        // IO Generator is done, shut down the IO workers...
                        backoff.spin();
                        backoff.spin(); // Very slight delay then issue terminate commands...

    //                        match rawseq_queue.push(ThreadCommand::Terminate) {
    //                            Ok(_) => (),
    //                            Err(x) => panic!("Unable to send command... {:#?}", x)
    //                        }
                        break;
                    }

                    match buffer[0] {
                        // 62 is a > meaning we have a new sequence id.
                        62 => {
                            jobs.fetch_add(1 as usize);
                            let wp = ThreadCommand::Work(seqbuffer[..seqlen].to_vec());
                            seqbuffer.clear();
                            seqlen = 0;

                            let mut result = rawseq_queue.push(wp);
                            while let Err(PushError(wp)) = result {
                                result = rawseq_queue.push(wp);
                            }

                            // pb.set_message(&format!("{}/2048 {}/131072", rawseq_queue.len(), seq_queue.len()));
                            pb.set_message(&msg);
                        },
                        _  => {
                            let slice_end = bytes_read.saturating_sub(1);
                            seqbuffer.extend_from_slice(&buffer[0..slice_end]);
                            seqlen = seqlen.saturating_add(slice_end);
                        }
                    }

                buffer.clear();
                }

            }
            *generator_done.write().unwrap() = true;
            for _ in 0..4 {
                let mut result = rawseq_queue.push(ThreadCommand::Terminate);
                while let Err(PushError(wp)) = result {
                    result = rawseq_queue.push(wp);
                }
            }}
            
            ).unwrap();
    }


    (seq_queue, jobs, generator_done, generator, children)
}

fn io_worker_thread(
    kmer_size: usize,
    rawseq_queue: Arc<ArrayQueue<ThreadCommand<Sequence>>>,
    seq_queue: Arc<ArrayQueue<ThreadCommand<Sequence>>>,
    jobs: Arc<AtomicCell<usize>>) {

    let backoff = Backoff::new();

    loop {
        if let Ok(command) = rawseq_queue.pop() {
            // We got the work packet!

            // We are finished, end the thread...
            if let ThreadCommand::Terminate = command {
                return;
            }

            let mut rawseq = command.unwrap();
            jobs.fetch_sub(1);

            capitalize_nucleotides(&mut rawseq);
            let coords = super::utils::get_good_sequence_coords(&rawseq);
            
            for (start_coords, end_coords) in coords {

                if (end_coords - start_coords) >= kmer_size {

                    jobs.fetch_add(1 as usize);
                    let wp = ThreadCommand::Work(rawseq[start_coords..end_coords].to_vec());

                    let mut result = seq_queue.push(wp);
                    while let Err(PushError(wp)) = result {
                        backoff.snooze();
                        result = seq_queue.push(wp);
                    }

                }
                
            }

        } else {
            backoff.snooze();
            backoff.reset();
        }
    }
}
