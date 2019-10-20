extern crate crossbeam;
extern crate thincollections;
extern crate wyhash;
extern crate once_cell;
extern crate opinionated;
extern crate bytecount;
extern crate flate2;
extern crate serde;
extern crate twox_hash;
// extern crate finalfrontier;
// pub extern crate finalfusion;
extern crate rand;
extern crate lru;

#[macro_use]
extern crate static_assertions;

mod utils;
pub mod kmer_counting;
mod threads;
pub mod vocab;
pub mod embeddings;
pub mod kmervocab;
