extern crate crossbeam;
extern crate thincollections;
extern crate wyhash;
extern crate once_cell;
extern crate bytecount;
extern crate flate2;
extern crate serde;
extern crate twox_hash;
// extern crate finalfrontier;
// pub extern crate finalfusion;
extern crate rand;
extern crate t1ha;
extern crate ndarray;
extern crate snap;
extern crate bincode;
extern crate lazy_static;
extern crate itertools;
extern crate linked_hash_set;
extern crate bumpalo;

// #[macro_use]
extern crate static_assertions; // I don't think this is being used anymore...

mod utils;
pub mod threads;
pub mod fasta;
pub mod kmers;
pub mod io;
pub mod sfasta;
pub mod gff3;
pub mod intervals;
