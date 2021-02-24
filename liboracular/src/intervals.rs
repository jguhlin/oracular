use hashbrown::HashMap;

use rust_lapper::Lapper;

// Landmark Interval Mapper

// Struct
// Landmark

pub struct IntervalMap<T: Eq + Clone + Sync + Send> {
    pub landmarks: HashMap<String, Lapper<u32, T>>,
}

impl<T: Eq + Clone + Sync + Send> Default for IntervalMap<T> {
    fn default() -> Self {
        IntervalMap {
            landmarks: Default::default(),
        }
    }
}

impl<T: Eq + Clone + Sync + Send> IntervalMap<T> {
    pub fn new() -> IntervalMap<T> {
        IntervalMap::default()
    }
}
