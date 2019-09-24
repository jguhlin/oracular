extern crate num_traits;
extern crate flate2;
extern crate indicatif;
extern crate opinionated;
extern crate seahash;
extern crate rayon;
extern crate twox_hash;
extern crate bitvec;
extern crate crossbeam;
extern crate dashmap;
extern crate fnv;
extern crate wyhash;
extern crate thincollections;
extern crate num_cpus;
extern crate rand;
extern crate rand_xorshift;

mod dictionary;

#[cfg(target_family = "unix")]
mod unix_specific {
    extern crate mimalloc;
    use mimalloc::MiMalloc;
    #[global_allocator]
    static GLOBAL: MiMalloc = MiMalloc;
}


use crossbeam::queue::{ArrayQueue, PushError};
use crossbeam::utils::Backoff;
use crossbeam::atomic::AtomicCell;

// Oracular is somehow a synonym of opinionated... a distant one, but still.
// That's where the name comes from.
// Opinionated is a minimizer based library
// Oracular is kmer / embedding based...

#[macro_use]
extern crate clap;
use clap::App;

use indicatif::ProgressBar;
use indicatif::ProgressStyle;

use std::sync::{Arc, RwLock};
use std::thread;
use std::thread::Builder;

use opinionated::fasta::{complement_nucleotides, capitalize_nucleotides};

use std::fs::File;
use std::io::prelude::*;
use std::io::{BufReader, Read};

/* Kmer size experiments...

capitalize before good seq
77.03user 4.61system 1:21.76elapsed 99%CPU (0avgtext+0avgdata 19030228maxresident)k

capitalize after good seq
76.83user 4.61system 1:21.55elapsed 99%CPU (0avgtext+0avgdata 19030136maxresident)k


without RC check...
76.88user 4.76system 1:21.75elapsed 99%CPU (0avgtext+0avgdata 19030152maxresident)k

with RC check...
76.83user 4.61system 1:21.55elapsed 99%CPU (0avgtext+0avgdata 19030136maxresident)k


WITH rc and proper capitalization

k=13
    342981308
     52270360
     6.561679

k=12
    343080666
     16353228
    20.979385

k=11
    343180024
      4235398
    81.026634

k=21
    342186444
    319348610
    1.0715138





// So the stats below (not timing, but kmer stats) are marred by an issue...

WITH reverse complement
   k=13 342981308
         66930112
         expected average: 5.12

   k=10 343279382 total
          5237227 unique
          expected average: 65.546

   k=11 343180024 total 
         11746932 unique
          expected average: 29.214

   k=9  343378740
          2235941
        expected average: 153.57


   These do not include reverse complement...
 * k=11 171590012 total tokens (This is what I suggest)
         10769668 unique kmers
         expected average: 15.93
   
   k=10 171639691
          4780515
         expected average: 35.9

   k=9  171689370 total
          2079595 unique
          expected avg: 82.559 (probably going towards useless?)

   k=25 171093222 total
        162723215 unique
        expected average: 1.05 (probably useless)

   k=17 171291938 total
        150872681 unique
        expected average: 1.135 (probably useless)

*/

/* Speed tests:
Just straight add
k=13
153.88user 11.19system 2:45.25elapsed 99%CPU (0avgtext+0avgdata 12136604maxresident)k

// With RwLock but single-threaded...
166.41user 12.44system 2:59.09elapsed 99%CPU (0avgtext+0avgdata 18407048maxresident)k
0inputs+0outputs (0major+10285501minor)pagefaults 0swaps

// With RwLock but multi-threaded...

sort and add based on previous id !!! ... so no
9683.21user 686.98system 5:49.42elapsed 2967%CPU (0avgtext+0avgdata 12135664maxresident)k
0inputs+0outputs (0major+126610995minor)pagefaults 0swaps

*/

const STACKSIZE: usize = 256 * 1024 * 1024;  // Stack size (needs to be > BUFSIZE + SEQBUFSIZE)
const WORKERSTACKSIZE: usize = 64 * 1024 * 1024;  // Stack size (needs to be > BUFSIZE + SEQBUFSIZE)

// type Sequences = Vec<Sequence>;
type Sequence = Vec<u8>;

enum ThreadCommand<T> {
    Work(T),
    Terminate,
}

impl ThreadCommand<Sequence> {
    // Consumes the ThreadCommand, which is just fine...
    fn unwrap(self) -> Sequence {
        match self {
            ThreadCommand::Work(x)   => x,
            ThreadCommand::Terminate => panic!("Unable to unwrap terminate command"),
        }
    }
}

fn main() {

    let yaml = load_yaml!("cli.yaml");
    let matches = App::from_yaml(yaml).get_matches();

    let kmer_size = value_t!(matches, "kmer", usize).unwrap_or(11);
    // let minn = value_t!(matches, "minn", usize).unwrap_or(13);
    // let maxn = value_t!(matches, "maxn", usize).unwrap_or(kmer_size.clone());
    // let step_size = value_t!(matches, "step", usize).unwrap_or(kmer_size.clone());
    // let w = value_t!(matches, "window", usize).unwrap_or(4);
    
    println!("k={}", kmer_size);

    // let test_file = "/mnt/data/nt/nt.gz";
    //  let test_file = "/mnt/data/3wasps/anno-refinement-run/genomes/Vvulg.fna";
    let test_file = "Vvulg.fna.gz";

    let jobs = Arc::new(AtomicCell::new(0 as usize));

/*
    let mut f = File::open(test_file).expect("Unable to open file");
    let mut byte_buffer = [0; 1024 * 1024 * 4];
    let mut buffer = Vec::new();
    
    let mut bytes_read: usize = 1;
    
    let now = Instant::now();
    let mut total: usize = 0;

    while bytes_read > 0 {
        bytes_read = f.read(&mut byte_buffer).expect("Unable to read file...");
        total += bytes_read;
    }
    println!("{} {}", now.elapsed().as_secs(), total);

    total = 0;

    drop(f);

    let mut f = File::open(test_file).expect("Unable to open file");
    buffer.push(IoSliceMut::new(&mut byte_buffer));

    bytes_read = 1;

    let now = Instant::now();

    while bytes_read > 0 {
        bytes_read = f.read_vectored(&mut buffer).expect("Unable to read file...");
        total += bytes_read;
    }
    println!("{} {}", now.elapsed().as_secs(), total);



    
    println!("Done"); */

/*    let pool = ThreadPool::new(32);
    let (tx, rx) = sync_channel(64);

    pool.execute(move|| {
        for mut kmers in rx.iter() {
            
        }        
    });

            pool.execute(move|| {
                        parse_kmers(tx, kmer_size, step_size, w, minn, maxn, slices_x);
                    });

    */

    // For testing at home, otherwise should be passable
    let num_threads = num_cpus::get();

    let seq_queue = Arc::new(ArrayQueue::<ThreadCommand<Sequence>>::new(1024 * 1));
    // let queue = Arc::new(ArrayQueue::<Vec<Vec<(usize, Vec<u8>)>>>::new(16));
    let done = Arc::new(RwLock::new(false));
    let generator_done = Arc::new(RwLock::new(false));

    // let aggregator;
    let generator;

    /* {
        
        let queue = Arc::clone(&queue);
        let done = Arc::clone(&done);
        

        aggregator = thread::Builder::new()
                            .name("Aggregator".to_string())
                            .spawn(move||
        {
            

            loop {
                // println!("Queue Length: {}", queue.len());
                while queue.is_empty() && !*done.read().unwrap() {
                    println!("");
                    println!("Waiting for data...");
                    println!("");
                    thread::sleep(time::Duration::from_millis(500));
                }

                if *done.read().unwrap() && queue.is_empty() {
                    println!("Finished, exiting... {} {:#?}", queue.len(), queue.is_empty());
                    break;
                }

                while !queue.is_empty() {
                    for result in queue.pop() {
                        for kmers in result {
                            for (hash, kmer) in kmers {
                                dict.add(hash, kmer);
                            }
                        }
                    }
                }
            }

            println!("{}", dict.tokens);
            println!("{}", dict.size);
            println!("{}", (dict.tokens as f32 / dict.size as f32));
        }).unwrap();

    } */

    let mut children = Vec::new();

    let dict_builder = match Builder::new()
                        .name("Dict Builder".into())
                        .spawn(|| Arc::new(dictionary::Dict::new()))
                    {
                        Ok(x)  => x,
                        Err(y) => panic!("{}", y)
                    };

    // let dict = Arc::new(dictionary::Dict::new());

    {
        // let queue = Arc::clone(&queue);
        let generator_done = Arc::clone(&generator_done);
        let seq_queue = Arc::clone(&seq_queue);
        let jobs = Arc::clone(&jobs);

        generator = thread::Builder::new()
                            .name("Generator".to_string())
                            .stack_size(STACKSIZE)
                            .spawn(move||
        {
            let filename = test_file.to_string();
            let mut buffer: Sequence = Vec::with_capacity(256);
            let mut seqbuffer: Sequence = Vec::with_capacity(2 * 1024 * 1024); // 2 Mb to start, will likely increase...
            let mut jobseqlen: usize = 0;
            let mut seqlen: usize = 0;
            // let pool = ThreadPool::new(48);

            let file = match File::open(&filename) {
                Err(why) => panic!("Couldn't open {}: {}", filename, why.to_string()),
                Ok(file) => file,
            };

            let pb = ProgressBar::new(file.metadata().unwrap().len());
            pb.set_style(ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {eta}")
                .progress_chars("█▇▆▅▄▃▂▁  "));

            // If file ends with .gz use flate2 to process it
            let fasta: Box<dyn Read> = if filename.ends_with("gz") {
                Box::new(flate2::read::GzDecoder::new(pb.wrap_read(file)))
            } else {
                Box::new(pb.wrap_read(file))
            };

            let mut reader = BufReader::with_capacity(32 * 1024 * 1024, fasta);
            let backoff = Backoff::new();

            while let Ok(bytes_read) = reader.read_until(b'\n', &mut buffer) {
                // File is empty, we are done!
                backoff.reset();
                if bytes_read == 0 {
                    break;
                }

                match buffer[0] {
                    // 62 is a > meaning we have a new sequence id.
                    // TODO: We don't care yet about the ID, but will soon...
                    62 => {
                        backoff.reset();
                        capitalize_nucleotides(&mut seqbuffer[..seqlen]);
                        // Rust's built-in uppercase function doesn't seem to save any time...
                        // seqbuffer[..seqlen].make_ascii_uppercase();
                        let coords = get_good_sequence_coords(&seqbuffer[..seqlen]);
                        
                        for (start_coords, end_coords) in coords {

                            // let mut push = Vec::new();
                            // push.push(seqbuffer[start_coords..end_coords].to_vec());
                            // jobseqlen.saturating_add(seqlen);

                            if (end_coords - start_coords) > (kmer_size * 2) {

                                jobs.fetch_add(1 as usize);
                                let wp = ThreadCommand::Work(seqbuffer[start_coords..end_coords].to_vec());

                                let mut result = seq_queue.push(wp);
                                while let Err(PushError(wp)) = result {
                                    backoff.snooze();
                                    result = seq_queue.push(wp);
                                }
                            }
                            
                        }
                        jobseqlen = jobseqlen.saturating_add(seqlen);
                        seqlen = 0;
                        seqbuffer.clear();

                        /* if seqlen > 0 {
                            work_packet.push((&seqbuffer[..seqlen]).to_vec());
                        }

                        jobseqlen = jobseqlen.saturating_add(seqlen);
                        seqlen = 0;
                        seqbuffer.clear();

                        if jobseqlen > JOBSIZE {
                            jobs.fetch_add(1);
                            let wp = ThreadCommand::Work(work_packet);

                            let mut result = seq_queue.push(wp);
                            while let Err(PushError(wp)) = result {
                                backoff.spin();
                                result = seq_queue.push(wp);
                            }

                            work_packet = Vec::new();
                            jobseqlen = 0;
                        } */

                    },
                    // Anything else is likely sequence we need...
                    _  => {
                        //seqbuffer.resize(seqlen + bytes_read - 1, 0);
                        // seqbuffer[seqlen..seqlen + bytes_read - 1].copy_from_slice(&buffer[0..bytes_read - 1]);
                        let slice_end = bytes_read.saturating_sub(1);
                        seqbuffer.extend_from_slice(&buffer[0..slice_end]);
                        seqlen = seqlen.saturating_add(slice_end);
                        // println!("{:#?}", String::from_utf8(seqbuffer[0..seqlen].to_vec()).unwrap());
                    }
                }

            buffer.clear();
            }

            *generator_done.write().unwrap() = true;
        }).unwrap();
    }

    let dict = dict_builder.join().unwrap();
    // let dict = OnceCell::new();

    for _ in 0..num_threads {
        let dict = Arc::clone(&dict);
        let seq_queue = Arc::clone(&seq_queue);
        let jobs = Arc::clone(&jobs);

        let child = match Builder::new()
                        .name("Worker".into())
                        .stack_size(WORKERSTACKSIZE)
                        .spawn(move || _worker_thread(kmer_size, dict, seq_queue, jobs)) {
                            Ok(x)  => x,
                            Err(y) => panic!("{}", y)
                        };
        
        children.push(child);
    }

    let backoff = Backoff::new();
    while !*generator_done.read().unwrap() {
        backoff.snooze();
    }

    println!("Generator done, joining...");

    generator.join().expect("Unable to join generator thread...");

    println!("Waiting for seq_queue to be empty...Currently at: {}", seq_queue.len());

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

    println!("Terminate commands sent, joining children");

    for child in children {
        match child.join() {
            Ok(_) => (),
            Err(x) => panic!("Error joining worker thread... {:#?}", x)
        }
    }

    println!("Children joined, getting dictionary stats...");

    *done.write().unwrap() = true;

    println!("{}", dict.tokens.load());
    println!("{}", dict.size.load());
    println!("{}", (dict.tokens.load() as f32 / dict.size.load() as f32));

}

fn get_good_sequence_coords (seq: &[u8]) -> Vec<(usize, usize)> {
    let mut start: Option<usize> = None;
    let mut end: usize;
    let mut cur: usize = 0;
    let mut start_coords;
    let mut end_coords;
    let mut coords: Vec<(usize, usize)> = Vec::new();
    let results = seq.windows(3).enumerate().filter(|(_y, x)| x != &[78, 78, 78]).map(|(y, _x)| y);
    for pos in results {
        match start {
            None    => { start = Some(pos); cur = pos; }
            Some(_x) => ()
        };

        if pos - cur > 1 {
            end = cur;
            start_coords = start.unwrap();
            end_coords = end;
            coords.push( (start_coords, end_coords) );
            start = None;
        } else {
            cur = pos;
        }
    }

    coords
}
fn _worker_thread(kmer_size: usize, 
                dict: Arc<dictionary::Dict>, 
                seq_queue: Arc<ArrayQueue<ThreadCommand<Sequence>>>,
                jobs: Arc<AtomicCell<usize>>) {

    loop {
        let backoff = Backoff::new();

        while seq_queue.is_empty() {
            backoff.snooze();
        }

        // Even when not empty, not guaranteed to be the thread to grab the work first...
        if let Ok(command) = seq_queue.pop() {
            backoff.reset();
            // We got the work packet!

            // We are finished, end the thread...
            if let ThreadCommand::Terminate = command {
                break;
            }          

            let rawseq = command.unwrap();
                jobs.fetch_sub(1);

                for i in 0..kmer_size {
                    rawseq[i..].chunks_exact(kmer_size).for_each(|x| { 
                        dict.add(&x);
                    });
                }

        }

    }
}
