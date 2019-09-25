extern crate crossbeam;
extern crate thincollections;
extern crate wyhash;
extern crate once_cell;
extern crate opinionated;
extern crate bytecount;

mod utils;
mod kmer_counting;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
