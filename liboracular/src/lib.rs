extern crate bytecount;
extern crate crossbeam;
extern crate flate2;
extern crate once_cell;
extern crate serde;
extern crate thincollections;
extern crate twox_hash;
extern crate wyhash;
// extern crate finalfrontier;
// pub extern crate finalfusion;
extern crate bincode;
extern crate bumpalo;
extern crate itertools;
extern crate lazy_static;
extern crate linked_hash_set;
extern crate ndarray;
extern crate rand;
extern crate snap;
extern crate t1ha;

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
