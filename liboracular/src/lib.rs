extern crate bumpalo;
extern crate bytecount;
extern crate crossbeam;
extern crate itertools;
extern crate jetscii;
extern crate lazy_static;
extern crate linked_hash_set;
extern crate ndarray;
extern crate once_cell;
extern crate rand;
extern crate snap;
extern crate thincollections;

// #[macro_use]
extern crate static_assertions; // I don't think this is being used anymore...

pub mod gff3;
pub mod intervals;
pub mod io;
pub mod kmers;
pub mod threads;
mod utils;
