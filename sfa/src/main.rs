extern crate mimalloc;

use liboracular::sfasta;
use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

extern crate indicatif;
extern crate rand;
extern crate rand_chacha;

use rand::prelude::*;
use rand_chacha::ChaCha20Rng;
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use indicatif::{HumanBytes, ProgressBar, ProgressIterator, ProgressStyle};

extern crate clap;
use clap::{load_yaml, App, ArgMatches};

fn main() {
    let yaml = load_yaml!("cli.yaml");
    let matches = App::from(yaml).get_matches();

    if let Some(matches) = matches.subcommand_matches("convert") {
        convert(&matches);
    }
    if let Some(matches) = matches.subcommand_matches("stats") {
        stats(&matches);
    }
    if let Some(matches) = matches.subcommand_matches("split") {
        split(&matches);
    }
    if let Some(matches) = matches.subcommand_matches("index") {
        index(&matches);
    }
}

fn convert(matches: &ArgMatches) {
    let fasta_filename = matches.value_of("input").unwrap();
    let path = Path::new(fasta_filename);
    sfasta::convert_fasta_file(
        &fasta_filename,
        &path.file_stem().unwrap().to_str().unwrap(),
    );
}

fn split(matches: &ArgMatches) {
    // TODO: Make the index as we go along, rather than at the end!
    let sfasta_filename = matches.value_of("input").unwrap();
    let output_prefix = matches.value_of("output").unwrap();
    let training_split = matches.value_of_t::<f32>("training").unwrap();
    let length_mode = matches.is_present("length_mode");
    let seed = matches.value_of_t::<u64>("seed").unwrap();

    let mut rng = ChaCha20Rng::seed_from_u64(seed);

    let out_train =
        File::create(output_prefix.to_string() + "_train.sfasta").expect("Unable to write to file");
    let mut out_train = BufWriter::new(out_train);

    let out_valid = File::create(output_prefix.to_string() + "_validation.sfasta")
        .expect("Unable to write to file");
    let mut out_valid = BufWriter::new(out_valid);

    let mut seqs = sfasta::Sequences::new(&sfasta_filename);
    // seqs.set_mode(sfasta::SeqMode::Random);
    let seqs = seqs.into_compressed_sequences();

    let mut total_len = 0;
    let mut training = 0;
    let mut validation = 0;
    // let mut validation_goal = 0;

    // Both get a copy of the header
    bincode::serialize_into(&mut out_train, &seqs.header)
        .expect("Unable to write to bincode output");
    bincode::serialize_into(&mut out_valid, &seqs.header)
        .expect("Unable to write to bincode output");

    // Split (based on number of sequences, or on total seq length)

    assert!(seqs.idx.is_some());
    let n = seqs.idx.as_ref().unwrap().0.len();

    let pb = ProgressBar::new(n as u64);

    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{bar:40.cyan/blue}] {eta_precise} {msg}")
            .progress_chars("█▛▌▖  ")
            .tick_chars("ACTGN"),
    );

    let mut seqs = seqs.progress_with(pb.clone());

    if !length_mode {
        // Can't do this without an index!

        let npct = (n as f32 * training_split) as usize;
        seqs.by_ref().take(npct).for_each(|s| {
            bincode::serialize_into(&mut out_train, &s).expect("Unable to write to bincode output")
        });
        seqs.take(n - npct).for_each(|s| {
            bincode::serialize_into(&mut out_valid, &s).expect("Unable to write to bincode output")
        });

        // TODO: I don't think this is complete / working...
        println!("Training: {}\tValidation: {}", npct, n - npct);
    } else {
        for entry in seqs {
            total_len += entry.len;
            let training_goal = (total_len as f32 * training_split) as usize;
            // validation_goal = (total_len as f32 * (1.0 - training_split)) as usize;
            let tchance = (training_goal as f32 - training as f32) / entry.len as f32;
            if rng.gen::<f32>() < tchance {
                bincode::serialize_into(&mut out_train, &entry)
                    .expect("Unable to write to bincode output");
                training += entry.len;
            } else {
                bincode::serialize_into(&mut out_valid, &entry)
                    .expect("Unable to write to bincode output");
                validation += entry.len;
            }

            pb.set_message(&format!(
                "Training Len: {} Validation Len: {}",
                HumanBytes(training),
                HumanBytes(validation)
            ));
        }

        println!(
            "Training Seq Length: {}\tValidation Seq Length:{}",
            training, validation
        );
    }

    // Drop (to close the files properly)
    drop(out_train);
    drop(out_valid);

    sfasta::index(&(output_prefix.to_string() + "_train.sfasta"));
    sfasta::index(&(output_prefix.to_string() + "_validation.sfasta"));
}

fn stats(matches: &ArgMatches) {
    let sfasta_filename = matches.value_of("input").unwrap();
    let (_, idx) = sfasta::open_file(&sfasta_filename);
    let mut seqs = sfasta::Sequences::new(&sfasta_filename);
    seqs.set_mode(sfasta::SeqMode::Random);
    let seqs = seqs.into_compressed_sequences();

    println!("Index Available: {}", idx.is_some());
    println!("Index Length: {}", idx.unwrap().0.len());
    println!("Header ID: {}", seqs.header.id.as_ref().unwrap());
    println!("Header Comment: {:#?}", seqs.header.comment);
    println!("Header Citation: {:#?}", seqs.header.citation);

    let ids: Vec<(String, u64)> = seqs.take(5).map(|seq| (seq.id, seq.len)).collect();
    println!("\nRandom (up to 5) IDs");
    for e in ids {
        println!("id: {} Length: {}", e.0, e.1);
    }

    let mut seqs = sfasta::Sequences::new(&sfasta_filename);
    seqs.set_mode(sfasta::SeqMode::Random);
    let mut seqs = seqs.into_compressed_sequences();
    let ec = seqs.next().unwrap();

    /*let seq = ec.clone().decompress();
    let len = if seq.len < 500 { seq.len as usize } else { 500 };
    println!("{}", std::str::from_utf8(&seq.seq[0..len]).unwrap());
    println!("{} {}", seq.len, ec.compressed_seq.len());*/
}

fn index(matches: &ArgMatches) {
    let sfasta_filename = matches.value_of("input").unwrap();
    sfasta::index(sfasta_filename);
}
