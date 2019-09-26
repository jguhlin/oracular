extern crate crossbeam;
extern crate thincollections;
extern crate wyhash;
extern crate once_cell;
extern crate opinionated;
extern crate bytecount;
extern crate bitpacking;
extern crate flate2;
extern crate serde;
extern crate bincode;
extern crate twox_hash;

mod utils;
pub mod kmer_counting;
mod threads;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
