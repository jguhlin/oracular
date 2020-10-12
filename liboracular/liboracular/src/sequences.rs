use std::fs::File;
use std::io::{BufReader, Read, BufRead};

use serde::{Serialize, Deserialize};

#[derive(PartialEq, Serialize, Deserialize)]
pub struct Sequence {
    pub seq: Vec<u8>,
    pub id:  String,
    pub taxons: Vec<usize>,
    pub taxon: usize,
}

