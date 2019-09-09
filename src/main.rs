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
extern crate finalfrontier;
extern crate finalfrontier_utils;
extern crate finalfusion;
extern crate finalfusion_utils;
extern crate rand;
extern crate rand_xorshift;

// Not convinced this actually works...
// Mimalloc compilation on windows is complicated. But
// the performance is worth it on unix.
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
// use std::ops::Not;

// use num_traits::{u64, u8, usize};

// use bitvec::prelude::*;

// Oracular is somehow a synonym of opinionated... a distant one, but still.
// That's where the name comes from.
// Opinionated is a minimizer based library
// Oracular is kmer / embedding based...

#[macro_use]
extern crate clap;
use clap::App;

use indicatif::ProgressBar;
use indicatif::ProgressStyle;

// use fnv::{FnvHashMap, FnvHashSet};

// type DnaKmerBinary = BitVec<BigEndian, u64>;

mod dictionary;
mod trainer;

//use std::fs::{OpenOptions};
//use std::io::{BufWriter};
// use std::hash::BuildHasherDefault;
// use std::collections::HashMap;
use std::sync::{Arc, RwLock};
// use std::sync::mpsc::{sync_channel};
use std::thread;
use std::thread::Builder;
// use std::time;

// use seahash::SeaHasher;

// use threadpool::ThreadPool;
// use rayon::prelude::*;
// use rayon::iter::ParallelBridge;
// use rayon::prelude::ParallelIterator;

// use opinionated::kmers::{Kmers, KmerOption};
use opinionated::fasta::{complement_nucleotides, capitalize_nucleotides};

use std::fs::File;
use std::io::prelude::*;
use std::io::{BufReader, Read};
// use std::time::{Duration, Instant};

/* fn mean(numbers: &[u32]) -> u32 {
    numbers.iter().sum::<u32>() as u32 / numbers.len() as u32
} */

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

// const BUFSIZE: usize = 128 * 1024 * 1024; // 128 Mb buffer size...
// const SEQBUFSIZE: usize = 8 * 1024 * 1024; // 8 Mb buffer for sequences
const STACKSIZE: usize = 256 * 1024 * 1024;  // Stack size (needs to be > BUFSIZE + SEQBUFSIZE)
const WORKERSTACKSIZE: usize = 64 * 1024 * 1024;  // Stack size (needs to be > BUFSIZE + SEQBUFSIZE)
//const JOBSIZE: usize = 512 * 1024; // 16 Mb, job size to send off to kmer processor.
                                         // Once buffered sequence length exceeds this, a job is sent
                                         // to the threadpool

// const JOBSIZE: usize = 2 * 1024 * 1024;

// type Sequences = Vec<Sequence>;
type Sequence = Vec<u8>;

enum ThreadCommand<T> {
    Work(T),
    Terminate,
}

/* impl ThreadCommand<Sequences> {
    // Consumes the ThreadCommand, which is just fine...
    fn unwrap(self) -> Sequences {
        match self {
            ThreadCommand::Work(x)   => x,
            ThreadCommand::Terminate => panic!("Unable to unwrap terminate command"),
        }
    }
} */

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

    // dict.set(dict_builder.join().unwrap());


/*
    let entries = opinionated::fasta::fasta_entries(test_file).unwrap();

    for i in entries {
        let mut x = match i {
                    Ok(x)  => x,
                    Err(y) => panic!("{}", y)
                };

        capitalize_nucleotides(&mut x.seq);
        let coords = get_good_sequence_coords(&x.seq);
        for (start_coords, end_coords) in coords {
            let seq = x.seq[start_coords..end_coords].to_vec();
            
            
            let kmers: Vec<Vec<u8>> = Kmers::with_step(&seq, kmer_size, 1)
                .filter(|x| 
                    match x {
                        KmerOption::Kmer(_) => true,
                        _                   => false
                })
                .map(|x| x.unwrap().to_vec()).into_iter().collect();

                let mut result = queue.push(kmers);
                while let Err(PushError(kmers)) = result {
                    let backoff = Backoff::new();
                    backoff.spin();
                    result = queue.push(kmers);
                }
*/
            /*
            kmers.par_sort();

            let mut previous_id: usize = 0;
            let mut previous_kmer: Vec<u8> = Vec::new();
            let mut previous_rc_id: usize = 0;

            for i in kmers {
                if i != previous_kmer {
                    previous_id = dict.add(&i);
                    
                    let mut rc = i.clone();
                    complement_nucleotides(&mut rc);
                    rc.reverse();

                    previous_rc_id = dict.add(&rc);
                    previous_kmer = i.to_vec();
                } else {
                    dict.add_to_id(previous_id);
                    dict.add_to_id(previous_rc_id);
                }

            } */



            
            // kmers.iter().for_each(|x| { dict.add(x) });
            /* kmers.iter().for_each(|x| { 
                let mut rc = x.to_vec();
                complement_nucleotides(&mut rc);
                rc.reverse();
                dict.add(&rc)
            }); */

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


    crate::trainer::create_embeddings(dict);

    
    // aggregator.join().expect("Unable to join aggregator thread...");
}




//    let pool = ThreadPool::new(32);

    // let (tx, rx) = sync_channel(64);

    /*let mut output = match matches.value_of("output") {
        None => { let output_fh = match OpenOptions::new().create(true).write(true).append(false).open("out.oracular") {
                Ok(x)  => x,
                Err(y) => panic!("{}", y)
            };
            BufWriter::new(output_fh) 
        },

        Some(filename) => { 
            let output_fh = match OpenOptions::new().create(true).write(true).append(false).open(filename) {
                Ok(x)  => x,
                Err(y) => panic!("{}", y)
            };

            BufWriter::new(output_fh)
        }};*/

    // let kmer_table: HashMap<Vec<u8>, Vec<u8>, BuildHasherDefault<SeaHasher>> = Default::default();

    // let mut total: usize = 0;
    // let mut total_jobs: u64 = 0;

    // let mut kmer_dict: FnvHashMap<Vec<u8>, u32> = Default::default();
    // let mut dict_kmer: Vec<Vec<u8>> = Vec::new();
    // let mut kmer_dict_id: u32 = 0;
    // let mut kmer_counts: FnvHashMap<Vec<u8>, usize> = Default::default();
    // let mut total_kmers: usize = 0;
    // let mut unique_kmers: FnvHashSet<Vec<u8>> = Default::default();

    /*

    {
        let entries = opinionated::fasta::fasta_entries(test_file).unwrap();

        for i in entries {
            let x = match i {
                        Ok(x)  => x,
                        Err(y) => panic!("{}", y)
                    };

            let coords = get_good_sequence_coords(&x.seq);
            let mut slices_x = Vec::new();
            let mut seqlength: usize = 0;
            for (start_coords, end_coords) in coords {
                slices_x.push(x.seq[start_coords..end_coords].to_vec());
                seqlength += end_coords - start_coords;

                if seqlength > 1024 * 1024 * 5 {
                    let tx = tx.clone();
                    pool.execute(move|| {
                        parse_kmers(tx, kmer_size, step_size, w, minn, maxn, slices_x);
                    });
                    total_jobs += 1;
                    slices_x = Vec::new();
                    seqlength = 0;
                }
            }

            if slices_x.len() > 0 {
                let tx = tx.clone();
                    pool.execute(move|| {
                        parse_kmers(tx, kmer_size, step_size, w, minn, maxn, slices_x);
                    });
                total_jobs += 1;
            }
        }
    } // Make sure entries gets dropped as soon as possible

    drop(tx); // Make sure the channel hangs up properly...

    // Single hash barely reduces memory and slows everything down

/*    let mut m: HashMap<u32, 
        HashMap<u32, usize, BuildHasherDefault<XxHash64>>,
        BuildHasherDefault<XxHash64>> = Default::default(); */
    let mut m: FnvHashMap<DnaKmerBinary, FnvHashMap<DnaKmerBinary, usize>> = Default::default();

    let pb = ProgressBar::new(total_jobs);
    pb.set_style(ProgressStyle::default_bar()
        .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {eta}")
        .progress_chars("█▇▆▅▄▃▂▁  "));


    // Entry API creates massive memory leak of over >100Gb
    // Or not.... hmm
    
    for mut hash in pb.wrap_iter(rx.iter()) {
        for (target, mut context_hash) in hash.drain() {
            let target_hash = m.entry(target).or_insert(Default::default());
            for (context,count) in context_hash.drain() {
                let ptr = target_hash.entry(context).or_insert(0);
                *ptr = *ptr + count;
                total += count;
            }
        // println!("Currently at: {} - {} - {}", total, pool.queued_count(), pool.active_count());
        }
    }

    println!("Finishing up");

    pool.join();

    */
// }

/*
fn parse_kmers( tx: std::sync::mpsc::SyncSender<FnvHashMap<DnaKmerBinary, FnvHashMap<DnaKmerBinary, usize>>>,
                kmer_size: usize, 
                step_size: usize, 
                w: usize,
                minn: usize,
                maxn: usize,
                xs: Vec<Vec<u8>>) 
        -> () {
        //-> Vec< (Arc<Vec<u8>>, Arc<Vec<u8>>) > {

    let mut collated: FnvHashMap<DnaKmerBinary, FnvHashMap<DnaKmerBinary, usize>> = Default::default();

    // println!("{}", xs.len());

    for mut x in xs {
        x.make_ascii_uppercase(); // I think this one is faster
        // TODO: Step size here is wrong for the context we are going for!
        // Need to change step size after processing entire sequence, not during...

        let kmers: Vec<_> = Kmers::with_step(&x, kmer_size, step_size)
            .filter(|x| 
                match x {
                    KmerOption::Kmer(_) => true,
                    _                   => false
                })
            .map(|x| x.unwrap().to_vec().clone()).into_iter().collect();

        let len = kmers.len();
        let mut start;
        let mut end;
        let mut target;
        let mut context: Vec<DnaKmerBinary> = Vec::new();
        let mut worked_on_counter: usize = 0;

        // println!("{}", len);

        for pos in 0..len {
            context.clear();
            worked_on_counter += 1;

            start = pos.saturating_sub(w);
            end   = pos.saturating_add(w);
            if end > len {
                end = len;
            }
                    
            target = convert_kmer_to_bits(&kmers[pos]);

            let mut target_hash = collated.entry(target).or_insert(Default::default());
            let mut ptr;

            if start < pos {
                for context_kmer in start..pos {
                    context.push(convert_kmer_to_bits(&kmers[context_kmer]));
                }
            }

            if end >= pos {
                for context_kmer in pos..end {
                    context.push(convert_kmer_to_bits(&kmers[context_kmer]));
                }
            }

            for ctx in &context {
                ptr = target_hash.entry((*ctx).clone()).or_insert(0);
                *ptr += 1;
                // data.push( (target.to_vec(), ctx.to_vec()) ); //(get_kmer(&kmer_table, target), get_kmer(&kmer_table, ctx) ) );
            }

            let target_kmer = (&kmers[pos]).to_vec();

            let mut rc = target_kmer.clone();
            complement_nucleotides(&mut rc);
            rc.reverse();

            // Errors are here...
            target_hash = collated.entry(convert_kmer_to_bits(&rc.to_vec())).or_insert(Default::default());

            for ctx in &context {
                ptr = target_hash.entry((*ctx).clone()).or_insert(0);
                *ptr += 1;
            }

            // Process subkmers

            for k in minn..=maxn {
                for subkmer in Kmers::with_step(&target_kmer, k, 1) {
                    match subkmer {
                        KmerOption::Kmer(subkmer) => {
                            target_hash = collated.entry(convert_kmer_to_bits(&subkmer)).or_insert(Default::default());
                            for ctx in &context {
                                ptr = target_hash.entry((*ctx).clone()).or_insert(0);
                                *ptr += 1;
                                //data.push( (subkmer.to_vec(), ctx.to_vec() ) ); //(get_kmer(&kmer_table, subkmer), get_kmer(&kmer_table, ctx) ) );
                            }
                        },
                        _ => () // Ignore partials, ignore empty
                    }
                }
            }

            for k in minn..=maxn { && x != &[110, 110, 110]
                for subkmer in Kmers::with_step(&rc, k, 1) {
                    match subkmer {
                        KmerOption::Kmer(subkmer) => {
                            target_hash = collated.entry(convert_kmer_to_bits(&subkmer)).or_insert(Default::default());
                            for ctx in &context {
                                ptr = target_hash.entry((*ctx).clone()).or_insert(0);
                                *ptr += 1;
                                //data.push( (subkmer.to_vec(), ctx.to_vec() )); //(get_kmer(&kmer_table, subkmer), get_kmer(&kmer_table, ctx) ) );
                            }
                        },
                        _ => () // Ignore partials, ignore empty
                    }
                }
            }

            // TODO: Do rc of contexts as well

/*            if worked_on_counter > 1024 * 512 { // 1024 * 1024 * 1 @ 50 is 28m
                                                // 1024 * 512 is 28m
                match tx.try_send(collated) {
                    Ok(_) => { println!("Sent..."); },
                    Err(x) => match x {
                        Full(x) => {
                            println!("Queue is full");
                            tx.send(x).expect("Unable to send hash-map to be collated");
                        }
                        _       => panic!("Ahhh!!!")
                    },
                    _  => panic!("Ahh!!!!!!")
                };
                // tx.send(collated).expect("Unable to send hash-map to be collated");
                collated = Default::default();
                worked_on_counter = 0;
            }
*/
    /*        if data.len() > 1024 * 512 {
                let collated = combine_entries(&mut data);
                tx.send(collated).expect("Unable to send data packet");
            } */
        }
    }
    
    if collated.len() > 0 {
        tx.send(collated).expect("Unable to send hash-map to be collated");
    }
}
*/
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

/*
#[inline(always)]
fn convert_seq_to_bits(seq: &[u8]) -> DnaKmerBinary {
    let mut binrep: DnaKmerBinary = BitVec::new();
    seq.iter().map(convert_to_bits)
        .for_each(|mut x| binrep.append(&mut x));
    binrep
        // .fold(BitVec::new(), |mut acc, mut x| { acc.append(&mut x); acc})
}

#[inline(always)]
fn bits_rc(seq: DnaKmerBinary) -> DnaKmerBinary {
    let mut seq = seq.not(); // Complement
    seq.reverse(); // And reverse
    seq
}

// Since we want to preserve N's, we can't use 2-bits. Trying out 3 bits here
#[inline(always)]
fn convert_to_bits(i: &u8) -> BitVec {
    match i {
        65 => bitvec![0,0,0], // A
        84 => bitvec![1,1,1], // T
        71 => bitvec![0,1,0], // G
        67 => bitvec![1,0,1], // C
        78 => bitvec![0,1,1], // N
        // [78] => bitvec![1,0,0], // N (for complement)
        // [78] => bitvec![1,1,0], // N (for reverse)
        // [78] => bitvec![0,0,1], // N (for reverse)

        97  => bitvec![0,0,0], // a -> A
        116 => bitvec![1,1,1], // t -> T
        103 => bitvec![0,1,0], // g -> G
        99  => bitvec![1,0,1], // c -> C
        110 => bitvec![0,1,1], // n -> N
        _        => panic!("Error, unknown nucleotide")
    }
} */


// 00000 AA !
// 00001 AT !
// 00010 GC !
// 00011 GG !
// 00100 AC !
// 00101 AG !
// 00110 CA !
// 00111 CT !
// 01000 AN !
// 01001 CN !
// 01010 NA !
// 01011 NC !
// Singletons for end of string
// 01100 A ! 
// 01101 C !
// 01110 N ! // Singleton
// 01111 NN // Placeholder, if "N" is flipped
// 10000 NN // Placeholder for double N's
// 10001 N ! // Singleton, if "N" is flipped
// 10010 G !
// 10011 T ! 
// 10100 NG !
// 10101 NT !
// 10110 GN !
// 10111 TN !
// 11000 GA !
// 11001 GT !
// 11010 TC !
// 11011 TG !
// 11100 CC !
// 11101 CG !
// 11110 TA !
// 11111 TT !

/*
#[inline(always)]
fn convert_to_bits_5bits(i: &[u8]) -> BitVec {
    match i {
        [65, 65] => bitvec![0,0,0,0,0], // AA complement is TT
        [84, 84] => bitvec![1,1,1,1,1], // TT -> AA
        [65, 84] => bitvec![0,0,0,0,1], // AT -> TA
        [84, 65] => bitvec![1,1,1,1,0], // TA -> AT
        [71, 67] => bitvec![0,0,0,1,0], // GC -> CG
        [67, 71] => bitvec![1,1,1,0,1], // CG -> GC
        [71, 71] => bitvec![0,0,0,1,1], // GG -> CC
        [67, 67] => bitvec![1,1,1,0,0], // CC -> GG
        [65, 67] => bitvec![0,0,1,0,0], // AC -> TG
        [84, 71] => bitvec![1,1,0,1,1], // TG -> AC
        [65, 71] => bitvec![0,0,1,0,1], // AG -> TC
        [84, 67] => bitvec![1,1,0,1,0], // TC -> AG
        [67, 65] => bitvec![0,0,1,1,0], // CA -> GT
        [71, 84] => bitvec![1,1,0,0,1], // GT -> CA
        [67, 84] => bitvec![0,0,1,1,1], // CT -> GA
        [71, 65] => bitvec![1,1,0,0,0], // GA -> CT
        [65, 78] => bitvec![0,1,0,0,0], // AN -> TN
        [84, 78] => bitvec![1,0,1,1,1], // TN -> AN
        [67, 78] => bitvec![0,1,0,0,1], // CN -> GN
        [71, 78] => bitvec![1,0,1,1,0], // GN -> CN
        [78, 65] => bitvec![0,1,0,1,0], // NA -> 
        [78, 84] => bitvec![1,0,1,0,1], // NT
        [78, 67] => bitvec![0,1,0,1,1], // NC
        [78, 71] => bitvec![1,0,1,0,0], // NG
        [78, 78] => bitvec![1,0,0,0,0], // NN
        // [78, 78] => bitvec![0,1,1,1,1], // NN
        [65]     => bitvec![0,1,1,0,0], // A
        [84]     => bitvec![1,0,0,1,1], // T
        [67]     => bitvec![0,1,1,0,1], // C
        [71]     => bitvec![1,0,0,1,0], // G
        [78]     => bitvec![1,0,0,0,1], // N
        // [78]     => bitvec![0,1,1,1,0], // N
        _        => panic!("Error, unknown doublet")
    }
} */

fn _worker_thread(kmer_size: usize, 
                dict: Arc<dictionary::Dict>, 
                seq_queue: Arc<ArrayQueue<ThreadCommand<Sequence>>>,
                jobs: Arc<AtomicCell<usize>>) {

    // let mut seq = Vec::new();
    // let mut results_packet;

    loop {
        //seq.clear();
        // results_packet = Vec::new();

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
            // for mut rawseq in command.unwrap() {
            // println!("Got job of size: {}", rawseq.len());
                jobs.fetch_sub(1);
                // capitalize_nucleotides(&mut rawseq);
                // let coords = get_good_sequence_coords(&rawseq);
                // for (start_coords, end_coords) in coords {
                    // seq.clear();
                    // seq.extend_from_slice(&rawseq[start_coords..end_coords]);
                    // seq.extend_from_slice(&rawseq[..]);
               
                // let kmers: Vec<Vec<u8>> = Kmers::with_step(&rawseq, kmer_size, 1)
                for i in 0..kmer_size {
                    rawseq[i..].chunks_exact(kmer_size).for_each(|x| { 
                        dict.add(&x);
                    });
                }

/*                let kmers: Vec<Vec<u8>> = Kmers::with_step(&rawseq, kmer_size, 1)
                    .filter(|x|
                        match x {
                        KmerOption::Kmer(_) => true,
                        _                   => false
                    })
                    .map(|x|
                    {
                        let kmer = x.unwrap().to_vec();
                        kmer
                    }).into_iter().collect();

                    for kmer in kmers {
                        dict.add(kmer);
                    } */
                //}
            // }


            // println!("Processed... {} {} remaining", jobs.load(), dict.tokens.load());
                    // println!("Added to dictionary, current total: {} {}", dict.size.load(), dict.tokens.load());
                    //results_packet.push(kmers);
                //}
            // }

            // let mut dict_w = dict.write().unwrap();

            //for kmers in results_packet {
                
            //}

            // println!("Added to dictionary, current total: {} {}", dict.size.load(), dict.tokens.load());
            // println!();
        }

    }
}
