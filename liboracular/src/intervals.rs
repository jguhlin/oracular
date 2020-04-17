use std::hash::BuildHasherDefault;
use std::collections::HashMap;
use twox_hash::XxHash64;

use rust_lapper::{Interval, Lapper};

// Landmark Interval Mapper

// Struct
// Landmark

pub struct IntervalMap<T: Eq + Clone> {
    pub landmarks: HashMap<String, 
                            Lapper<T>, 
                            BuildHasherDefault<XxHash64>>,
}

impl<T: Eq + Clone> Default for IntervalMap<T> {
    fn default() -> Self { 
        IntervalMap {
            landmarks: Default::default(),
        }
    }
}

impl<T: Eq + Clone> IntervalMap<T> {
    pub fn new() -> IntervalMap<T> {
        IntervalMap::default()
    }
}

