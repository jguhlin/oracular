extern crate bincode;
extern crate bumpalo;
extern crate bytecount;
extern crate crossbeam;
extern crate flate2;
extern crate itertools;
extern crate lazy_static;
extern crate linked_hash_set;
extern crate ndarray;
extern crate once_cell;
extern crate rand;
extern crate rayon;
extern crate serde;
extern crate serde_bytes;
extern crate snap;
extern crate t1ha;
extern crate thincollections;
extern crate twox_hash;
extern crate wyhash;
extern crate zstd;

// #[macro_use]
extern crate static_assertions; // I don't think this is being used anymore...

pub mod fasta;
pub mod gff3;
pub mod intervals;
pub mod io;
pub mod kmers;
pub mod sfasta;
pub mod threads;
mod utils;
