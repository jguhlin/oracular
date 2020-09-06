extern crate num_traits;
extern crate flate2;
extern crate indicatif;
// extern crate rayon;
extern crate crossbeam;
extern crate fnv;
extern crate wyhash;
extern crate thincollections;
extern crate num_cpus;
extern crate liboracular;
extern crate serde;
extern crate snap;
extern crate bincode;

// extern crate mimalloc;
// use mimalloc::MiMalloc;
// #[global_allocator]
// static GLOBAL: MiMalloc = MiMalloc;

// Oracular is (very distantly) a synonym of opinionated...
// That's where the name comes from.
// Opinionated is a minimizer based library
// Oracular is kmer / embedding based...

#[macro_use]
extern crate clap;
use clap::App;

fn main() {

    let yaml = load_yaml!("cli.yaml");
    let _matches = App::from_yaml(yaml).get_matches();

/*     if let Some(matches) = matches.subcommand_matches("generate-embeddings") {
        generate_embeddings(matches);
    } else
    if let Some(matches) = matches.subcommand_matches("query") {
        query_embeddings(matches);
    } */

    println!("I don't do anything right now");
}

