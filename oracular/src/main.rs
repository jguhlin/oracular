extern crate num_traits;
extern crate flate2;
extern crate indicatif;
extern crate opinionated;
// extern crate rayon;
extern crate crossbeam;
extern crate fnv;
extern crate wyhash;
extern crate thincollections;
extern crate num_cpus;
extern crate liboracular;
extern crate snap;

extern crate mimalloc;
use mimalloc::MiMalloc;
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

// Oracular is (very distantly) a synonym of opinionated...
// That's where the name comes from.
// Opinionated is a minimizer based library
// Oracular is kmer / embedding based...

#[macro_use]
extern crate clap;
use clap::App;

use std::fs::File;
use std::io::{BufReader, Read, BufRead};
use liboracular::vocab::build_vocab_from_finaldict;

use finalfrontier::{
    SentenceIterator, SimpleVocab, SkipgramTrainer, SubwordVocab, Vocab, VocabBuilder,
    WriteModelBinary, SGD, SubwordVocabConfig, NGramConfig
};

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
    let test_file = "/mnt/data/3wasps/anno-refinement-run/genomes/Vvulg.fna";
    // let test_file = "Vvulg.fna.gz";

    let num_threads = num_cpus::get();
    let filename = test_file.clone();

    /* let pb = ProgressBar::new(file.metadata().unwrap().len());
    pb.set_style(ProgressStyle::default_bar()
        .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {eta}")
        .progress_chars("█▇▆▅▄▃▂▁  "));

    let file = match File::open(&filename) {
        Err(why) => panic!("Couldn't open {}: {}", filename, why.to_string()),
        Ok(file) => file,
    }; */

    // If file ends with .gz use flate2 to process it

    println!("Counting kmers with {} threads", num_threads);
    let dict = liboracular::kmer_counting::count_kmers(num_threads, kmer_size, filename);

    println!("{}", dict.tokens.load());
    println!("{}", dict.size.load());
    println!("{}", (dict.tokens.load() as f32 / dict.size.load() as f32));

    let final_dict = dict.convert_to_final();

    let app = SkipGramApp::new();

    let vocab: SubwordVocab<_, _> = build_vocab_from_finaldict(final_dict);

    // let mut out_fh = snap::Writer::new(File::create(format!("{}.bc", "vvulg")).unwrap());
    // bincode::serialize_into(&mut out_fh, &final_dict).expect("Unable to write to bincode file");
}