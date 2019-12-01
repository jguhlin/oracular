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
extern crate serde;
extern crate snap;
extern crate bincode;


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

use std::path::Path;
use std::fs::File;
use std::io::BufReader;


use liboracular::vocab::build_vocab_from_finaldict;
use liboracular::embeddings::train;
use liboracular::model::Model;


fn main() {

    let yaml = load_yaml!("cli.yaml");
    let matches = App::from_yaml(yaml).get_matches();

    if let Some(matches) = matches.subcommand_matches("generate-embeddings") {
        generate_embeddings(matches);
    } else
    if let Some(matches) = matches.subcommand_matches("query") {
        query_embeddings(matches);
    }
}

fn query_embeddings(matches: &clap::ArgMatches<'_>) {
    let kmer_size = value_t!(matches, "kmer", usize).unwrap_or(11);
    let model_filename = value_t!(matches, "input", String).expect("Invalid input identified");
    let vocab_filename = value_t!(matches, "vocab", String).expect("Invalid input identified");

    let model = Model::new(model_filename, vocab_filename, kmer_size);

    println!("Model loaded, querying...");

    let seq = "TAAGACTGTTCGGAAGTACTAGAGATATGGTTTCTTTGCGAATCTAAGTTATCGTGGTTAATGTCCTTCTCATATCGCTTTTTGTATGACTTGTTGCAAGTAGGACGACGTACCTAAGATATTTTTATTTATTTTAACTTATATTAGGATCACATCGATGTGAGAGTATTAATTGGAAAAGAAAGAAAGATAAAAGGAACACGTTTCATCCACTCTACATATAAATATTATTGCGATATTCACCTTCCGTTTTTCTTGGATGTTTTTATGAAAATCGTGCGACGGGCTTTGTAGGAACGATAATGTAGGTTTATTTCTATGCTTGGAAACTTGCATCATCGTAGGTTTAACATTTGGAAAAGAGCAAATCCTTGGAGACATTTTTGAGGATATGGATTTGAAAGTAGGAAGCATTTTATTTTGTGTACACAACATTATATCTGGTTTAATTTCGTGTCCCACTTTTTGAAGCTCACTTTTAGCGTACCTTTTTTTTCTTAACATTTTGTTACTCATTTCGATATTTTGAAGTTGTATAAAATGTAGCTAGGAATAGGATAAAAGAACAATTATTCAACTTTCGACTCGTACAACTGATACGTTATACATAATATTGCGATTCTATTTAACGCAAAATCATTCTGATACGATAATTACTTGTACTCACCAAAGTTAACAGACTTGAAACATAAACGTAGTTACGATACGTGAAGATTGCGGTTTTTGGCGTAGCAGCTTTAGATAAGACTATTGAATGAAATATTTTATTAAAACTTTATAAATAAAATGTGAGAAACGATTATTTTGCCATGAAAAATCAAGTTACTCCTACGATTAAACAAAAATGTTCAAACGATATGGTTTCGTTAAGTTTCGTGCAGCCTTGTGAATTAGCCTATTATAAATATATAAAACTTTATGATACGTACGTATACGTACTTACGAGACATTACTTACCTGAAAGTTTCTTGGGTATCTACTACATCACGGAAATTTGGTTAAACATTTTTTCGAGAGAAAAAAATGTCATAGGCCACATGAGATTTTACTATAAGTGAGCCTTCGTGAGCCAGAGGTCGTACTGGATCGTATGTACGCATATACAAGTATACCCATTAGTTTCAAGTTCCATGGTTGTAAGATTTTTCTCACGTATTGAATAGATAAACATAGAAAAGAGTCACGACGGTGGGAAAAACGAATGGAATAATCGATGTCAATTAAACAAATATCGGAGAAAACGTGAAAGAATACAAATTTTATGGTAATAGTTTGAAAATAAAAATTATATTTCTATATACGCATATAAGTAGGAAAGAAAATATATCGAAGTAATATTATGTGCAAGTATTTACATATATAGCAAATACCTATGAGATCAAAGAATATTCTTCTGCAAATGATAAACGTATATTATCATTGTTTTTCCGGTAACCATATGTCAGATGTGCAATCACCACCAGGGTATCTGTTTCGAGAAGAATTGGATTATTGATCGATCGATCGATCTTAATTGGGTGTATCGTAAATTaaaagtaacgcgaaaaaac";
    let result0 = model.get_weighted_embedding(seq);

    let seq = "GAATAGGATAAAAGAACAATTATTCAACTTTCGACTCGTACAACTGATACGTTATACATAATATTGCGATTCTATTTAACGCAAAATCATTCTGATACGATAATTACTTGTACTCACCAAAGTTAACAGACTTGAAACATAAACGTAGTTACGATACGTGAAGATTGCGGTTTTTGGCGTAGCAGCTTTAGATAAGACTA";
    let result1 = model.get_weighted_embedding(seq);

    let seq = "AGGCGACCACACTGGCGCGGGAAATGGGTTATACCGAACCGGACCCGCGAGATGATCTTTCTGGTATGGAAGGCGACCACACTGGCGCGGGAAATGGGTTATACCGAACCGGACCCGCGAGATGATCTTTCTGGTATGGAAGGCGACCACACTGGCGCGGGAAATGGGTTATACCGAACCGGACCCGCGAGATGATCTTTCTGGTATGGAAGGCGACCACACTGGCGCGGGAAATGGGTTATACCGAACCGGACCCGCGAGATGATCTTTCTGGTATGGAAGGCGACCACACTGGCGCGGGAAATGGGTTATACCGAACCGGACCCGCGAGATGATCTTTCTGGTATGGAAGGCGACCACACTGGCGCGGGAAATGGGTTATACCGAACCGGACCCGCGAGATGATCTTTCTGGTATGGAAGGCGACCACACTGGCGCGGGAAATGGGTTATACCGAACCGGACCCGCGAGATGATCTTTCTGGTATGGA";
    let result2 = model.get_weighted_embedding(seq);

    let seq = "TTTATTATAATTATGATTACTTATCACGACGCATTCGCGAAAGCGAACAATTACCTTGATGATGCAAATCTTTATTATAATTATGATTACTTATCACGACGCATTCGCGAAAGCGAACAATTACCTTGATGATGCAAATCTTTATTATAATTATGATTACTTATCACGACGCATTCGCGAAAGCGAACAATTACCTTGATGATGCAAATCTTTATTATAATTATGATTACTTATCACGACGCATTCGCGAAAGCGAACAATTACCTTGATGATGCAAATCTTTATTATAATTATGATTACTTATCACGACGCATTCGCGAAAGCGAACAATTACCTTGATGATGCAAATCTTTATTATAATTATGATTACTTATCACGACGCATTCGCGAAAGCGAACAATTACCTTGATGATGCAAATCTTTATTATAATTATGATTACTTATCACGACGCATTCGCGAAAGCGAACAATTACCTTGATGATGCAAATC";
    let result3 = model.get_weighted_embedding(seq);

    let cos_result = Model::cos_distance(&result0, &result1);
    println!("Cos: {:#?}", cos_result);

    let cos_result = Model::cos_distance(&result0, &result2);
    println!("Cos: {:#?}", cos_result);

    let cos_result = Model::cos_distance(&result2, &result3);
    println!("Cos: {:#?}", cos_result);

    println!("Result0: {:#?}", result0);
    println!("Result1: {:#?}", result1);
    println!("Result2: {:#?}", result2);
    println!("Result3: {:#?}", result3);


/*
    let embeddings0 = model.get_weighted_embedding("CCTAAGCTACT");
    let embeddings1 = model.get_weighted_embedding("CCTAAGCTACT");
    let cos_result = Model::cos_distance(&embeddings0, &embeddings1);
    println!("Should be 1 Cos: {:#?}", cos_result); */

/*
    let embeddings0 = model.get_embeddings_and_count(b"CCTAAGCTACT");
    let embeddings1 = model.get_embeddings_and_count(b"AGCTGATCGAG");

    let result = &embeddings0.0.dot(&embeddings1.0);
    println!("Cos: {:#?}", result);

    let embeddings2 = model.get_embeddings_and_count(b"TGTAAGGACCG");

    // let result = (embeddings0.0 - embeddings1.0).dot(&embeddings1.0).sqrt();
    let result = euclidian_distance(&embeddings0.0, &embeddings1.0);
    println!("{:#?}", result);

    // let result = embeddings0.0.dot(&embeddings2.0).sqrt();
    let result = euclidian_distance(&embeddings0.0, &embeddings2.0);
    println!("{:#?}", result);

    let result = model.embeddings.embedding_similarity(model.embeddings.embedding("CCTAAGCTACT").unwrap().as_view(), 5);
    println!("{:#?}", result); */

}


fn generate_embeddings(matches: &clap::ArgMatches<'_>) {
    let kmer_size = value_t!(matches, "kmer", usize).unwrap_or(13);
    // let minn = value_t!(matches, "minn", usize).unwrap_or(13);
    // let maxn = value_t!(matches, "maxn", usize).unwrap_or(kmer_size.clone());
    // let step_size = value_t!(matches, "step", usize).unwrap_or(kmer_size.clone());
    let context_size = value_t!(matches, "window", usize).unwrap_or(5);
    let num_threads = value_t!(matches, "threads", usize).unwrap_or(16);
    let dims = value_t!(matches, "dims", usize).unwrap_or(32);
    let epochs = value_t!(matches, "epochs", usize).unwrap_or(5);
    let mincount = value_t!(matches, "mincount", usize).unwrap_or(5);
    let min_n = value_t!(matches, "min_n", usize).unwrap_or(9);
    let max_n = value_t!(matches, "max_n", usize).unwrap_or(11);
    
    println!("k={}", kmer_size);

    // let test_file = "/mnt/data/nt/nt.gz";
    // let test_file = "/mnt/data/3wasps/anno-refinement-run/genomes/Vvulg.fna";
    // let test_file = "Vvulg.fna.gz";

    let filename = value_t!(matches, "input", String).expect("Invalid input identified");

    // let num_threads = num_cpus::get();
    // let filename = test_file.clone();

    /* let pb = ProgressBar::new(file.metadata().unwrap().len());
    pb.set_style(ProgressStyle::default_bar()
        .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {eta}")
        .progress_chars("█▇▆▅▄▃▂▁  "));

    let file = match File::open(&filename) {
        Err(why) => panic!("Couldn't open {}: {}", filename, why.to_string()),
        Ok(file) => file,
    }; */

    // If file ends with .gz use flate2 to process it

    let corpus_path = Path::new(&filename);
    let filename_stem = corpus_path.file_stem().unwrap().to_str().unwrap();
    let output_filename = format!("kmer_counts_{}_k{}.bc", filename_stem, kmer_size.to_string());

    let final_dict;

    if Path::new(&output_filename).exists() {
        println!("Existing kmer count file found -- opening...");
        let kmercounts_fh = snap::Reader::new(BufReader::new(File::open(output_filename.clone()).unwrap()));
        final_dict = bincode::deserialize_from(&mut BufReader::new(kmercounts_fh)).expect("Unable to read to bincode file");
        println!("Finished reading file");
    } else {
        println!("Counting kmers with {} threads", num_threads);
        let dict = liboracular::kmer_counting::count_kmers(num_threads, kmer_size, &filename);
        final_dict = dict.convert_to_final();
        let mut kmercounts_fh = snap::Writer::new(File::create(output_filename.clone()).unwrap());
        bincode::serialize_into(&mut kmercounts_fh, &final_dict).expect("Unable to write to bincode file");
        println!("Saved kmer count results to {}", output_filename);
    }

    println!("Dictionary - Total Tokens: {}", final_dict.tokens);
    println!("Dictionary - Total Kmers: {}", final_dict.entries);
    println!("Dictionary - Avg. occurences of kmers: {}", (final_dict.tokens as f32 / final_dict.entries as f32));
    println!();

    // let app = SkipGramApp::new();

    let vocab = build_vocab_from_finaldict(final_dict, mincount, min_n, max_n);
    let model = train(num_threads, dims, epochs, context_size, vocab, &filename, kmer_size);

    // let mut out_fh = snap::Writer::new(File::create(format!("{}.bc", "vvulg")).unwrap());
    // bincode::serialize_into(&mut out_fh, &final_dict).expect("Unable to write to bincode file");
}
