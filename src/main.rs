extern crate flate2;
extern crate indicatif;
extern crate opinionated;
extern crate seahash;
extern crate rayon;
extern crate twox_hash;
extern crate fnv;
extern crate bitvec;

use twox_hash::XxHash64;
use std::sync::mpsc::TrySendError::Full;


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

    let (tx, rx) = sync_channel(64);

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

    let mut total: usize = 0;
    let mut total_jobs: u64 = 0;

    // let mut kmer_dict: FnvHashMap<Vec<u8>, u32> = Default::default();
    // let mut dict_kmer: Vec<Vec<u8>> = Vec::new();
    // let mut kmer_dict_id: u32 = 0;
    // let mut kmer_counts: FnvHashMap<Vec<u8>, usize> = Default::default();
    // let mut total_kmers: usize = 0;
    // let mut unique_kmers: FnvHashSet<Vec<u8>> = Default::default();

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
}

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

            for k in minn..=maxn {
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

#[inline(always)]
fn convert_kmer_to_bits(kmer: &[u8]) -> DnaKmerBinary {
    let mut binrep: DnaKmerBinary = BitVec::new();
    kmer.chunks(2)
        .map(convert_to_bits)
        .for_each(|mut x| binrep.append(&mut x));
    binrep
        // .fold(BitVec::new(), |mut acc, mut x| { acc.append(&mut x); acc})
}

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
        // [78, 78] => bitvec![0,1,1,1,1], // NN
        [65]     => bitvec![0,1,1,0,0], // A
        [84]     => bitvec![1,0,0,1,1], // T
        [67]     => bitvec![0,1,1,0,1], // C
        [71]     => bitvec![1,0,0,1,0], // G
        [78]     => bitvec![1,0,0,0,1], // N
        // [78]     => bitvec![0,1,1,1,0], // N
        _        => panic!("Error, unknown doublet")
    }
}