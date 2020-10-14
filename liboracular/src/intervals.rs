use hashbrown::HashMap;

use rust_lapper::Lapper;

// Landmark Interval Mapper

// Struct
// Landmark

pub struct IntervalMap<T: Eq + Clone> {
    pub landmarks: HashMap<String, Lapper<T>>,
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
