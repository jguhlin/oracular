extern crate flate2;
extern crate indicatif;
extern crate opinionated;
extern crate seahash;
extern crate rayon;
extern crate twox_hash;
extern crate fnv;
extern crate bitvec;

use twox_hash::XxHash64;

#[macro_use]
use bitvec::prelude::*;

// Oracular is somehow a synonym of opinionated... a distant one, but still.
// That's where the name comes from.
// Opinionated is a minimizer based library
// Oracular is kmer / embedding based...

#[macro_use]
extern crate clap;
use clap::App;

use indicatif::ProgressBar;
use indicatif::ProgressStyle;

use fnv::{FnvHashMap, FnvHashSet};

type DnaKmerBinary = BitVec<BigEndian, u64>;

//use std::fs::{OpenOptions};
//use std::io::{BufWriter};
use std::hash::BuildHasherDefault;
use std::collections::HashMap;
use std::sync::{Arc, RwLock};
use std::sync::mpsc::{sync_channel};

// use seahash::SeaHasher;

use threadpool::ThreadPool;
use rayon::prelude::*;
// use rayon::iter::ParallelBridge;
// use rayon::prelude::ParallelIterator;

use opinionated::kmers::{Kmers, KmerOption};
use opinionated::fasta::{complement_nucleotides};

fn mean(numbers: &[u32]) -> u32 {
    numbers.iter().sum::<u32>() as u32 / numbers.len() as u32
}

fn main() {
    let yaml = load_yaml!("cli.yaml");
    let matches = App::from_yaml(yaml).get_matches();

    let kmer_size = value_t!(matches, "kmer", usize).unwrap_or(15);
    let minn = value_t!(matches, "minn", usize).unwrap_or(13);
    let maxn = value_t!(matches, "maxn", usize).unwrap_or(kmer_size.clone());
    let step_size = value_t!(matches, "step", usize).unwrap_or(kmer_size.clone());
    let w = value_t!(matches, "window", usize).unwrap_or(4);
    
    let test_file = "/mnt/data/3wasps/anno-refinement-run/genomes/Vvulg.fna";

    let pool = ThreadPool::new(32);

    let (tx, rx) = sync_channel(16);

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

    let kmer = "ATGTTGGCAAACAAGAAATACTAAGTGGTCTGGA";
    let kmer2 = "ATGTTGGCAAACAAGAAATACTAAGTGGTCTGGAA";
    let kmer3 = "ATGTTGGCAAACAAGAAATACTAAGTGGTCTGGAAA";

    println!("{:#?}", convert_kmer_to_bits(kmer));
    println!("{:#?}", convert_kmer_to_bits(kmer2));
    println!("{:#?}", convert_kmer_to_bits(kmer3));

    let mut total: usize = 0;
    let mut total_jobs: u64 = 0;

    let mut kmer_dict: FnvHashMap<Vec<u8>, u32> = Default::default();
    let mut dict_kmer: Vec<Vec<u8>> = Vec::new();
    let mut kmer_dict_id: u32 = 0;
    let mut kmer_counts: FnvHashMap<Vec<u8>, usize> = Default::default();
    let mut total_kmers: usize = 0;
    let mut unique_kmers: FnvHashSet<Vec<u8>> = Default::default();

    {
        let entries = opinionated::fasta::fasta_entries(test_file).unwrap();
        for i in entries {
            println!("Processed");
            let x = match i {
                Ok(x)    => x,
                Err(y)  => panic!("{}", y)
            };

            let result: (FnvHashSet<Vec<u8>>, usize) = 
                get_good_sequence_coords(&x.seq)
                .par_iter()
                .map(|(start_coords, end_coords)| {
                    let mut slice = x.seq[*start_coords..*end_coords].to_vec();
                    slice.make_ascii_uppercase();

                    let kmers = get_all_kmers(&mut slice, &kmer_size, &minn, &maxn, &step_size);
                    let mut unique_kmers: FnvHashSet<Vec<u8>> = Default::default();
                    let mut total_count: usize = kmers.len();
                    unique_kmers.extend(kmers);
                    (unique_kmers, total_count)
                })
                .reduce(
                    || (Default::default(), 0 as usize),
                    |mut a, b| {
                            (a.0).extend(b.0); 
                            (a.0, a.1 + b.1)
                    });

            unique_kmers.extend(result.0);
            total_kmers += result.1;

        }

        dict_kmer.resize(unique_kmers.len(), Vec::new()); // Probably speeds things up

        unique_kmers
            .iter()
            .enumerate()
            .for_each(|(id, kmer)|
                {
                    kmer_dict.insert(kmer.clone(), id as u32);
                    dict_kmer[id] = kmer.clone()
                });
    }

    println!("{}", kmer_dict_id);

    println!("{}", String::from_utf8(kmer_dict.keys().next().unwrap().to_vec()).unwrap());

    drop(unique_kmers); // Don't need this anymore...
    drop(kmer_counts); // Not using this one just yet...



    // let keep = counts.clone().into_iter().enumerate().filter(|&(_x, y)| y >= 5).collect::<Vec<_>>();
    // counts = counts.into_iter().filter(|&x| x >= 5).collect::<Vec<u32>>();

    // for (x,y) in keep {
    //    if y > 10000 {
    //        let kmer = String::from_utf8(dict_kmer[x as usize].clone()).unwrap();
    //        println!("{} {}", kmer, y);
    //    }
    //}
    
    // println!("Most occurences of a single kmer: {}", counts.iter().max().unwrap());
    // println!("Fewest occurences of a single kmer: {}", counts.iter().min().unwrap());
    // println!("Mean occurences of kmers: {}", mean(&counts));

    println!("Found {} unique kmers", kmer_dict.len());
    println!("Found {} total kmers", total_kmers);
//    println!("Total after removing those with fewer than 5: {}", counts.len());
   
    let kmer_dict = Arc::new(kmer_dict);
    let dict_kmer = Arc::new(dict_kmer);

    {
        let entries = opinionated::fasta::fasta_entries(test_file).unwrap();

        for i in entries {
            // data.clear();
            let x = match i {
                        Ok(x)  => x,
                        Err(y) => panic!("{}", y)
                    };

            let coords = get_good_sequence_coords(&x.seq);
            for (start_coords, end_coords) in coords {
                let slice = x.seq[start_coords..end_coords].to_vec();

                let kmer_dict = Arc::clone(&kmer_dict);
                let dict_kmer = Arc::clone(&dict_kmer);

                let tx = tx.clone();
                pool.execute(move|| {
                    parse_kmers(tx, kmer_size, step_size, w, minn, maxn, slice, kmer_dict, dict_kmer);
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
    let mut m: FnvHashMap<u32, FnvHashMap<u32, usize>> = Default::default();

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
}

fn parse_kmers( tx: std::sync::mpsc::SyncSender<FnvHashMap<u32, FnvHashMap<u32, usize>>>,
                kmer_size: usize, 
                step_size: usize, 
                w: usize,
                minn: usize,
                maxn: usize,
                mut x: Vec<u8>,
                kmer_dict: Arc<FnvHashMap<Vec<u8>,u32>>,
                dict_kmer: Arc<Vec<Vec<u8>>>) 
        -> () {
        //-> Vec< (Arc<Vec<u8>>, Arc<Vec<u8>>) > {

    let mut collated: FnvHashMap<u32, FnvHashMap<u32, usize>> = Default::default();

    x.make_ascii_uppercase(); // I think this one is faster
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
    let mut context: Vec<u32> = Vec::new();
    // let mut data: Vec< (Vec<u8>, Vec<u8>) > = Vec::new();
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
                
        target = kmer_dict.get(&kmers[pos]).expect("Kmer is missing! 1");

        let mut target_hash = collated.entry(*target).or_insert(Default::default());
        let mut ptr;

        if start < pos {
            for context_kmer in start..pos {
                context.push(*kmer_dict.get(&kmers[context_kmer]).expect("Kmer is missing! 2"));
            }
        }

        if end >= pos {
            for context_kmer in pos..end {
                context.push(*kmer_dict.get(&kmers[context_kmer]).expect("Kmer is missing! 3"));
            }
        }

        for ctx in &context {
            ptr = target_hash.entry(*ctx).or_insert(0);
            *ptr += 1;
            // data.push( (target.to_vec(), ctx.to_vec()) ); //(get_kmer(&kmer_table, target), get_kmer(&kmer_table, ctx) ) );
        }

        let target_kmer = (&kmers[pos]).to_vec();

        let mut rc = target_kmer.clone();
        complement_nucleotides(&mut rc);
        rc.reverse();

        // Errors are here...
        target_hash = collated.entry(*kmer_dict.get(&rc.to_vec()).expect("Kmer is missing! 4")).or_insert(Default::default());

        for ctx in &context {
            ptr = target_hash.entry(ctx.clone()).or_insert(0);
            *ptr += 1;
        }

        // Process subkmers

        for k in minn..=maxn {
            for subkmer in Kmers::with_step(&target_kmer, k, 1) {
                match subkmer {
                    KmerOption::Kmer(subkmer) => {
                        target_hash = collated.entry(*kmer_dict.get(subkmer).expect("Kmer is missing! 5")).or_insert(Default::default());
                        for ctx in &context {
                            ptr = target_hash.entry(*ctx).or_insert(0);
                            *ptr += 1;
                            //data.push( (subkmer.to_vec(), ctx.to_vec() ) ); //(get_kmer(&kmer_table, subkmer), get_kmer(&kmer_table, ctx) ) );
                        }
                    },
                    _ => () // Ignore partials, ignore empty
                }
            }
        }

        for k in minn..=maxn {
            for subkmer in Kmers::with_step(&rc, k, 1) {
                match subkmer {
                    KmerOption::Kmer(subkmer) => {
                        target_hash = collated.entry(*kmer_dict.get(subkmer).expect("Kmer is missing! 6")).or_insert(Default::default());
                        for ctx in &context {
                            ptr = target_hash.entry(*ctx).or_insert(0);
                            *ptr += 1;
                            //data.push( (subkmer.to_vec(), ctx.to_vec() )); //(get_kmer(&kmer_table, subkmer), get_kmer(&kmer_table, ctx) ) );
                        }
                    },
                    _ => () // Ignore partials, ignore empty
                }
            }
        }

        // TODO: Do rc of contexts as well

        if worked_on_counter > 1024 * 512 {
            tx.send(collated).expect("Unable to send hash-map to be collated");
            collated = Default::default();
            worked_on_counter = 0;
        }

/*        if data.len() > 1024 * 512 {
            let collated = combine_entries(&mut data);
            tx.send(collated).expect("Unable to send data packet");
        } */
    }
    
    if collated.len() > 0 {
        tx.send(collated).expect("Unable to send hash-map to be collated");
    }
}

fn get_good_sequence_coords (seq: &[u8]) -> Vec<(usize, usize)> {
    let mut start: Option<usize> = None;
    let mut end: usize;
    let mut cur: usize = 0;
    let mut start_coords;
    let mut end_coords;

    let mut coords: Vec<(usize, usize)> = Vec::new();

    let results = seq.windows(3).enumerate().filter(|(_y, x)| (x != &[78, 78, 78] && x != &[110, 110, 110])).map(|(y, _x)| y);
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

fn get_all_kmers(seq: &mut [u8], &kmer_size: &usize, &minn: &usize, &maxn: &usize, &step_size: &usize) -> Vec<Vec<u8>> {
    let mut result: Vec<Vec<u8>> = Vec::new();

    seq.make_ascii_uppercase();

    // Get all kmers in reverse complement of sequence as well
    let mut rc = seq.to_vec();
    complement_nucleotides(&mut rc);
    rc.reverse();

    let kmers = Kmers::with_step(&seq, kmer_size, step_size)
                    .chain(Kmers::with_step(&rc, kmer_size, step_size));

    // Get subkmers
    let maxnmer = if maxn == kmer_size {
            maxn - 1
        } else {
            maxn
        };

    for k in minn..=maxnmer {
        // Step size has to be 1 for subkmers
        for kmer in Kmers::with_step(&seq, k, 1) {
            if let KmerOption::Kmer(k) = kmer {
                result.push(k.to_vec());
            }
        }

        for kmer in Kmers::with_step(&rc,  k, 1)  {
            if let KmerOption::Kmer(k) = kmer {
                result.push(k.to_vec());
            }
        }
    }

    for kmer in kmers {
        if let KmerOption::Kmer(k) = kmer {
            result.push(k.to_vec());
        }
    }

    result
}

#[inline(always)]
fn convert_kmer_to_bits(kmer: &str) -> DnaKmerBinary {
    kmer.as_bytes()
        .chunks(2)
        .map(convert_to_bits)
        .fold(BitVec::new(), |mut acc, mut x| { acc.append(&mut x); acc})
}


#[inline(always)]
fn convert_to_bits(i: &[u8]) -> BitVec {
    match i {
        [65, 65] => bitvec![0,0,0,0,0], // AA
        [84, 84] => bitvec![1,1,1,1,1], // TT
        [65, 84] => bitvec![0,0,0,0,1], // AT
        [84, 65] => bitvec![1,1,1,1,0], // TA
        [71, 67] => bitvec![0,0,0,1,0], // GC
        [67, 71] => bitvec![1,1,1,0,1], // CG
        [71, 71] => bitvec![0,0,0,1,1], // GG
        [67, 67] => bitvec![1,1,1,0,0], // CC
        [65, 67] => bitvec![0,0,1,0,0], // AC
        [84, 71] => bitvec![1,1,0,1,1], // TG
        [65, 71] => bitvec![0,0,1,0,1], // AG
        [84, 67] => bitvec![1,1,0,1,0], // TC
        [67, 65] => bitvec![0,0,1,1,0], // CA
        [71, 84] => bitvec![1,1,0,0,1], // GT
        [67, 84] => bitvec![0,0,1,1,1], // CT
        [71, 65] => bitvec![1,1,0,0,0], // GA
        [65, 78] => bitvec![0,1,0,0,0], // AN
        [84, 78] => bitvec![1,0,1,1,1], // TN
        [67, 78] => bitvec![0,1,0,0,1], // CN
        [71, 78] => bitvec![1,0,1,1,0], // GN
        [78, 65] => bitvec![0,1,0,1,0], // NA
        [78, 84] => bitvec![1,0,1,0,1], // NT
        [78, 67] => bitvec![0,1,0,1,1], // NC
        [78, 71] => bitvec![1,0,1,0,0], // NG
        [78, 78] => bitvec![1,0,0,0,0], // NN
        [65]     => bitvec![0,1,1,0,0], // A
        [84]     => bitvec![1,0,0,1,1], // T
        [67]     => bitvec![0,1,1,0,1], // C
        [71]     => bitvec![1,0,0,1,0], // G
        _        => panic!("Error, unknown doublet")
    }
}