extern crate crossbeam;
extern crate mimalloc;
extern crate rand;
extern crate rand_xoshiro;

use crossbeam::queue::ArrayQueue;
use crossbeam::utils::Backoff;
use mimalloc::MiMalloc;
use std::sync::atomic::Ordering;
use std::thread;
use std::thread::{park, JoinHandle};

use rand::prelude::*;
use std::sync::atomic::AtomicBool;
use std::sync::Arc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use rand_xoshiro::rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

// NOTE: New naming convention
// Rust-y stuff is "iter" Python is "Generator"

use liboracular::kmers::KmerWindow;
use liboracular::kmers::{
    rc_kmerwindow, replace_random, DiscriminatorMasked, DiscriminatorMaskedGenerator, Gff3Kmers,
    Gff3KmersIter, KmerCoordsWindow, KmerCoordsWindowIter, KmerWindowGenerator,
};
use liboracular::sfasta;
use liboracular::io as io;

use pyo3::prelude::*;
// use pyo3::wrap_pyfunction;
use pyo3::types::PyDict;
use pyo3::PyIterProtocol;

// use pyo3::wrap_pyfunction;

use pyo3::wrap_pyfunction;

// use std::fs::File;
// use std::io::BufReader;

// TODO: Should be literals to concatenate...
#[inline]
fn convert_string_to_array(k: usize, s: &[u8]) -> Vec<bool> {
    let mut out: Vec<bool> = vec![false; k * 5];

    for (x, c) in s.iter().enumerate() {
        match c {
            65 => out[5*x] = true,   // A
            84 => out[5*x+1] = true, // T
            67 => out[5*x+3] = true, // C
            71 => out[5*x+4] = true, // G
            78 => out[5*x+2] = true, // N
            _  => out[5*x+2] = true, // N for everything else...
            // A_  => { out[5*x+2] = 1; println!("Invalid Character! {} in {}", c, std::str::from_utf8(s).unwrap()) }  // N
        };
    }

    out
}

// ** Kmer Classification GFF3
#[pyclass]
struct Gff3KmerGenerator {
    iter: Box<dyn Iterator<Item = Gff3Kmers> + Send>,
    k: usize,
    offset: usize,
    window_size: usize,
    filename: String,
    gff3filename: String,
    rc: bool,
    types: Vec<String>,
    rand: bool,
}

#[pyproto]
impl PyIterProtocol for Gff3KmerGenerator {
    fn __iter__(mypyself: PyRefMut<Self>) -> PyResult<PyObject> {
        let gil = Python::acquire_gil();
        let py = gil.python();
        Ok(mypyself.into_py(py))
    }

    fn __next__(mut mypyself: PyRefMut<Self>) -> PyResult<Option<PyObject>> {
        let mut finished = false;
        let mut item = None;

        while !finished {
            item = match mypyself.iter.next() {
                Some(x) => {
                    finished = true;
                    Some(x)
                }
                None => {
                    mypyself.offset += 1;
                    if mypyself.k == mypyself.offset && mypyself.rc {
                        println!("Finished, at the correct step...");
                        return Ok(None);
                    } else {
                        if mypyself.k == mypyself.offset {
                            mypyself.rc = true;
                            mypyself.offset = 0;
                        }

                        let kmercoords_window_iter = KmerCoordsWindowIter::new(
                            &mypyself.filename,
                            mypyself.k,
                            mypyself.window_size,
                            mypyself.offset,
                            mypyself.rc,
                            mypyself.rand,
                        );

                        let iter = Gff3KmersIter::new(
                            &mypyself.gff3filename,
                            kmercoords_window_iter,
                            mypyself.k,
                        );

                        mypyself.iter = Box::new(iter);
                        continue;
                    }
                }
            };
        }

        match item {
            Some(x) => {
                let Gff3Kmers {
                    kmers,
                    classifications,
                    id,
                    coords,
                    rc,
                } = x;
                let kmers: Vec<Vec<bool>> = kmers
                    .iter()
                    .map(|x| convert_string_to_array(mypyself.k, x))
                    .collect();

                let gil = Python::acquire_gil();
                let py = gil.python();
                let pyout = PyDict::new(py);
                // let pyout = PyTuple::new(py, [kmers, truth]);
                pyout.set_item("kmers", kmers).expect("Py Error");
                pyout
                    .set_item("classifications", classifications)
                    .expect("Py Error");
                pyout.set_item("id", id).expect("Py Error");
                pyout.set_item("rc", rc).expect("Py Error");
                pyout.set_item("coords", coords).expect("Py Error");
                Ok(Some(pyout.to_object(py)))
            }
            None => Ok(None),
        }
    }
}

#[pymethods]
impl Gff3KmerGenerator {
    #[new]
    fn new(
        k: usize,
        filename: String,
        window_size: usize,
        gff3filename: String,
        rand: bool,
    ) -> Self {
        // Create KmerWindowGenerator
        let kmercoords_window_iter =
            KmerCoordsWindowIter::new(&filename, k, window_size, 0, false, rand);

        let iter = Gff3KmersIter::new(&gff3filename, kmercoords_window_iter, k);
        let types = iter.types.clone();

        Gff3KmerGenerator {
            iter: Box::new(iter),
            k,
            offset: 0,
            filename,
            gff3filename,
            window_size,
            types,
            rc: false,
            rand,
        }
    }

    fn types(&self) -> Vec<String> {
        self.types.clone()
    }
}

// ** Queue-version of DiscriminatorMaskedGeneratorWrapper
// Get around the PYTHON GIL
#[pyclass]
struct MaskedKmersGenerator {
    queueimpl: QueueImpl<WindowSubmission>,
}

type Kmer = Vec<bool>;
type Id = String;
type Truth = bool;

type WindowKmers = Vec<Kmer>;
type WindowIds = Id;
type WindowTruths = Vec<Truth>;

type WindowSubmission = (WindowKmers, WindowIds, WindowTruths);

#[pymethods]
impl MaskedKmersGenerator {
    #[new]
    fn new(
        k: usize,
        filename: String,
        window_size: usize,
        replacement_pct: f32,
        rand: bool,
        queue_size: usize,
    ) -> Self {
        let queueimpl = QueueImpl::new(queue_size, 16, move |shutdown, exhausted, queue| {
            let mut offset = 0;
            let mut rc = false;

            loop {
                // Create KmerWindowGenerator
                let kmer_window_generator =
                    KmerWindowGenerator::new(&filename, k, window_size, offset, rc, rand);

                let discriminator_masked_generator =
                    DiscriminatorMaskedGenerator::new(replacement_pct, k, kmer_window_generator);

                let mut iter = Box::new(discriminator_masked_generator);

                while let Some(x) = iter.next() {
                    if shutdown.load(Ordering::Relaxed) {
                        return; // We are done, something triggered a
                                // shutdown...
                    }

                    let DiscriminatorMasked { kmers, id, truth } = x;
                    let kmers = kmers
                        .iter()
                        .map(|x| convert_string_to_array(k, x))
                        .collect();

                    let mut batch = (kmers, id, truth);
                    while let Err(x) = queue.push(batch) {
                        // Test if we are prematurely shutdown...
                        if shutdown.load(Ordering::Relaxed) {
                            return; // We are done, something triggered a
                                    // shutdown...
                        }
                        batch = x;
                        park(); // Queue is full, park the thread...
                    }
                }

                offset += 1;

                if (k == offset) && rc {
                    return;
                } else if k == offset {
                    offset = 0;
                    rc = true;
                }
            }
        });

        MaskedKmersGenerator { queueimpl }
    }

    fn len(mypyself: PyRef<Self>) -> usize {
        mypyself.queueimpl.queue.len()
    }
}

#[pyproto]
impl PyIterProtocol for MaskedKmersGenerator {
    fn __iter__(mypyself: PyRefMut<Self>) -> PyResult<PyObject> {
        let gil = Python::acquire_gil();
        let py = gil.python();
        Ok(mypyself.into_py(py))
    }

    fn __next__(mut mypyself: PyRefMut<Self>) -> PyResult<Option<PyObject>> {
        if mypyself.queueimpl.is_finished() {
            return Ok(None);
        }

        // Unpark the thread...
        mypyself.queueimpl.unpark();

        let mut result = mypyself.queueimpl.queue.pop();
        let backoff = Backoff::new();

        while result == None {
            mypyself.queueimpl.unpark();
            backoff.snooze();

            // Check for exhaustion (or shutdown)...
            if mypyself.queueimpl.is_finished() {
                return Ok(None);
            }

            result = mypyself.queueimpl.queue.pop();
        }

        let result = result.unwrap();
        let gil = Python::acquire_gil();
        let py = gil.python();
        let pyout = PyDict::new(py);
        pyout
            .set_item("kmers", result.0)
            .expect("Error with Python");
        pyout.set_item("id", result.1).expect("Error with Python");
        pyout
            .set_item("truth", result.2)
            .expect("Error with Python");

        // One last unparking...
        mypyself.queueimpl.unpark();

        Ok(Some(pyout.to_object(py)))
    }
}

/// ** MatchedKmersGenerator
/// 1 if both sets of sequences are from the same sequence, 0 if not

#[pyclass]
struct MatchedKmersGenerator {
    queueimpl: QueueImpl<MatchedSubmission>,
}

type MatchedKmers = (WindowKmers, WindowKmers);
type Matches = bool;

type MatchedSubmission = (MatchedKmers, Matches);

#[pymethods]
impl MatchedKmersGenerator {
    #[new]
    fn new(k: usize, filename: String, window_size: usize, queue_size: usize) -> Self {
        let queueimpl = QueueImpl::new(queue_size, 16, move |shutdown, exhausted, queue| {
            let mut offset = 0;
            let mut rc = false;
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

            // TODO: Make even smarter -- Load up 1k windows and pick from there matching
            // and non-matching ones, including some RC ones as well...
            loop {
                // Create KmerWindowGenerator
                let mut iter1 =
                    KmerWindowGenerator::new(&filename, k, window_size, offset, rc, true);

                let mut iter2 =
                    KmerWindowGenerator::new(&filename, k, window_size, offset, rc, true);

                loop {
                    if shutdown.load(Ordering::Relaxed) {
                        return; // We are done, something triggered a
                                // shutdown...
                    }

                    let mut item1;
                    let mut item2;
                    let matched;

                    if rng.gen::<bool>() {
                        matched = true;

                        item1 = match iter1.next() {
                            Some(x) => x,
                            None => break,
                        };

                        // Half the time skip a window...
                        if rng.gen::<bool>() {
                            iter1.next();
                        }

                        item2 = match iter1.next() {
                            Some(x) => x,
                            None => break,
                        };

                        while item1.id != item2.id {
                            item1 = item2.clone();

                            item2 = match iter1.next() {
                                Some(x) => x,
                                None => break,
                            };
                        }
                    } else {
                        matched = false;

                        item1 = match iter1.next() {
                            Some(x) => x,
                            None => break,
                        };

                        item2 = match iter1.next() {
                            Some(x) => x,
                            None => break,
                        };

                        while item1.id == item2.id {
                            item2 = match iter2.next() {
                                Some(x) => x,
                                None => break,
                            };
                        }
                    }

                    let KmerWindow {
                        kmers: kmers1,
                        id: _,
                        rc: _,
                    } = item1;

                    let KmerWindow {
                        kmers: kmers2,
                        id: _,
                        rc: _,
                    } = item2;

                    let kmers1 = kmers1
                        .iter()
                        .map(|x| convert_string_to_array(k, x))
                        .collect();

                    let kmers2 = kmers2
                        .iter()
                        .map(|x| convert_string_to_array(k, x))
                        .collect();

                    let mut batch = ((kmers1, kmers2), matched);
                    while let Err(x) = queue.push(batch) {
                        // Test if we are prematurely shutdown...
                        if shutdown.load(Ordering::Relaxed) {
                            return; // We are done, something triggered a
                                    // shutdown...
                        }
                        batch = x;
                        park(); // Queue is full, park the thread...
                    }
                }

                offset += 1;

                if (k == offset) && rc {
                    return;
                } else if k == offset {
                    offset = 0;
                    rc = true;
                }
            }
        });

        MatchedKmersGenerator { queueimpl }
    }

    fn len(mypyself: PyRef<Self>) -> usize {
        mypyself.queueimpl.queue.len()
    }
}

#[pyproto]
impl PyIterProtocol for MatchedKmersGenerator {
    fn __iter__(mypyself: PyRefMut<Self>) -> PyResult<PyObject> {
        let gil = Python::acquire_gil();
        let py = gil.python();
        Ok(mypyself.into_py(py))
    }

    fn __next__(mut mypyself: PyRefMut<Self>) -> PyResult<Option<PyObject>> {
        if mypyself.queueimpl.is_finished() {
            return Ok(None);
        }

        // Unpark the thread...
        mypyself.queueimpl.unpark();

        let mut result = mypyself.queueimpl.queue.pop();
        let backoff = Backoff::new();

        while result == None {
            mypyself.queueimpl.unpark();
            backoff.snooze();

            // Check for exhaustion (or shutdown)...
            if mypyself.queueimpl.is_finished() {
                return Ok(None);
            }

            result = mypyself.queueimpl.queue.pop();
        }

        let result = result.unwrap();
        let gil = Python::acquire_gil();
        let py = gil.python();
        let pyout = PyDict::new(py);
        pyout
            .set_item("kmers", result.0)
            .expect("Error with Python");
        pyout
            .set_item("matched", result.1)
            .expect("Error with Python");

        // One last unparking...
        mypyself.queueimpl.unpark();

        Ok(Some(pyout.to_object(py)))
    }
}

/// ** TripleLossKmersGenerator
/// 3 outputs
///   [] -> is kmer replaced (up to 15% replaced, or whatever is passed in),
///   1 if both sets of sequences are from the same sequence, 0 if not
///   0 if not reverse complement, 1 if reverse complement of first sequence
/// (thus the second output is going to be 1 as well!) All can have up to 15%
/// replaced

#[pyclass]
struct TripleLossKmersGenerator {
    queueimpl: QueueImpl<TripleLossSubmission>,
}

type ReverseComplement = bool;
type Truths = Vec<bool>;

type TripleLossSubmission = (MatchedKmers, (Truths, Truths, Matches, ReverseComplement));

#[pymethods]
impl TripleLossKmersGenerator {
    #[new]
    fn new(
        k: usize,
        filename: String,
        replacement_pct: f32,
        window_size: usize,
        queue_size: usize,
    ) -> Self {
        let queueimpl = QueueImpl::new(queue_size, 16, move |shutdown, exhausted, queue| {
            let mut offset = 0;
            let mut rc = false;
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

            let mut sfasta = sfasta::Sequences::new(&filename);

            // let headers = sfasta.idx.as_ref().unwrap().0.clone();
            let locs = sfasta.idx.as_ref().unwrap().1.clone();

            // TODO: Make even smarter -- Load up 1k windows and pick from there matching
            // and non-matching ones, including some RC ones as well...

            loop {
                let mut count = 0;
                // Create KmerWindowGenerator
                let mut iter1 =
                    KmerWindowGenerator::new(&filename, k, window_size, offset, rc, true);

                'inner: loop {
                    if shutdown.load(Ordering::Relaxed) {
                        return; // We are done, something triggered a
                                // shutdown...
                    }

                    let mut item1;
                    let mut item2;
                    let matched;
                    let reversecomplement;

                    // So we have 3 states..
                    // Matched sequence -- So kmer window 1 and 2 from the same sequence
                    // Matched Sequence -- RC -- Kmer window 1 and it's reverse complement
                    // Non-matched sequence -- Kmer window 1 and 2 from completely different seqs...
                    let choice: u8 = rng.gen_range(0, 3); // Give us a number between 0 and 2

                    // Always need a starting window...
                    item1 = match iter1.next() {
                        Some(x) => x,
                        None => {
                            break 'inner;
                        }
                    };

                    count += 1;

                    while item1.kmers.len() < window_size {
                        item1 = match iter1.next() {
                            Some(x) => x,
                            None => {
                                break 'inner;
                            }
                        };
                    }

                    // Matched Sequence
                    if choice == 0 {
                        matched = true;
                        reversecomplement = false;

                        item2 = match get_random_sequence_from_id(
                            &mut sfasta,
                            k,
                            window_size,
                            &item1.id,
                            &mut rng,
                        ) {
                            Some(x) => x,
                            None => continue,
                        };

                    // RC
                    } else if choice == 1 {
                        matched = true;
                        reversecomplement = true;

                        item2 = item1.clone();
                        item2 = rc_kmerwindow(item2);

                    // Not matching sequence...
                    } else {
                        matched = false;
                        reversecomplement = false;

                        item2 = match get_random_sequence_from_locs(
                            &mut sfasta,
                            k,
                            window_size,
                            &locs,
                            &mut rng,
                        ) {
                            Some(x) => x,
                            None => continue,
                        };
                    }

                    let KmerWindow {
                        kmers: mut kmers1,
                        id: _,
                        rc: _,
                    } = item1;
                    let KmerWindow {
                        kmers: mut kmers2,
                        id: _,
                        rc: _,
                    } = item2;

                    let truth1 = replace_random(k, replacement_pct, &mut kmers1, &mut rng);
                    let truth2 = replace_random(k, replacement_pct, &mut kmers2, &mut rng);

                    let kmers1 = kmers1
                        .iter()
                        .map(|x| convert_string_to_array(k, x))
                        .collect();

                    let kmers2 = kmers2
                        .iter()
                        .map(|x| convert_string_to_array(k, x))
                        .collect();

                    let mut batch = (
                        (kmers1, kmers2),
                        (truth1, truth2, matched, reversecomplement),
                    );

                    while let Err(x) = queue.push(batch) {
                        // Test if we are prematurely shutdown...
                        if shutdown.load(Ordering::Relaxed) {
                            return; // We are done, something triggered a
                                    // shutdown...
                        }
                        batch = x;
                        park(); // Queue is full, park the thread...
                    }
                }

                println!("{} {} Total from Iter: {}", offset, rc, count);

                offset += 1;

                if (k == offset) && rc {
                    return;
                } else if k == offset {
                    offset = 0;
                    rc = true;
                }

                println!("Offset: {} RC: {}", offset, rc);
            }
        });

        TripleLossKmersGenerator { queueimpl }
    }

    fn len(mypyself: PyRef<Self>) -> usize {
        mypyself.queueimpl.queue.len()
    }
}

/// Support functions for triple loss generator
fn is_all_ns(seq: &[u8]) -> bool {
    if bytecount::count(&seq, b'N') == seq.len() {
        true
    } else {
        false
    }
}

/// Support functions for triple loss generator
fn get_random_sequence_from_id<R: Rng + ?Sized>(
    sfasta: &mut sfasta::Sequences,
    k: usize,
    window_size: usize,
    id: &str,
    rng: &mut R,
) -> Option<KmerWindow> {
    let needed_length = (k * window_size) + k;

    let mut seq;

    seq = sfasta.get(&id).unwrap();
    if seq.seq.len() < needed_length {
        return None;
    }

    let seqlen = seq.seq.len().saturating_sub(needed_length);

    if seqlen == 0 {
        return None;
    }

    let mut start = rng.gen_range(0, seqlen);
    let mut end = start + needed_length;
    assert!(end < seq.seq.len());

    while is_all_ns(&seq.seq[start..end]) {
        start = rng.gen_range(0, seqlen);
        end = start + needed_length;
    }

    seq.seq = seq.seq[start..end].to_vec();

    let mut iter2 = match KmerWindowGenerator::from_sequence(
        seq,
        k,
        window_size,
        0, // Already a random sequence, random offset won't do anything...
        rng.gen(),
    ) {
        Some(x) => x,
        None => return None,
    };

    iter2.next()
}

/// Support functions for triple loss generator
fn get_random_sequence_from_locs<R: Rng + ?Sized>(
    sfasta: &mut sfasta::Sequences,
    k: usize,
    window_size: usize,
    locs: &Vec<u64>,
    mut rng: &mut R,
) -> Option<KmerWindow> {
    let needed_length = ((k * window_size) + k) as u64;

    let mut seqlen = 0;

    let mut loc = locs.choose(&mut rng).unwrap().clone();
    let mut seq = sfasta.get_header_at(loc).expect("Unable to get header");

    while (seq.len < needed_length) {
        loc = locs.choose(&mut rng).unwrap().clone();
        seq = sfasta.get_header_at(loc).unwrap();
    }

    let seqlen = seq.len.saturating_sub(needed_length);

    let mut start = rng.gen_range(0, seqlen);
    let mut end = start + needed_length;
    assert!(end < seq.len);
    if seq.len >= end + 1000 {
        end += 1000;
    } else {
        end = seq.len;
    }

    let sequence = sfasta.get_seq_slice(seq.id.clone(), start, end).expect("Unable to get seq slice");

    let mut workseq = io::Sequence {
        seq: sequence,
        end: end as usize,
        location: start as usize,
        id: seq.id.clone()
    };

    let mut iter2 = match KmerWindowGenerator::from_sequence(
        workseq,
        k,
        window_size,
        0, // Already a random sequence, random offset won't do anything...
        rng.gen(),
    ) {
        Some(x) => x,
        None => return None,
    };

    iter2.next()
}

#[pyproto]
impl PyIterProtocol for TripleLossKmersGenerator {
    fn __iter__(mypyself: PyRefMut<Self>) -> PyResult<PyObject> {
        let gil = Python::acquire_gil();
        let py = gil.python();
        Ok(mypyself.into_py(py))
    }

    fn __next__(mut mypyself: PyRefMut<Self>) -> PyResult<Option<PyObject>> {
        if mypyself.queueimpl.is_finished() {
            return Ok(None);
        }

        // Unpark the thread...
        mypyself.queueimpl.unpark();

        let mut result = mypyself.queueimpl.queue.pop();
        let backoff = Backoff::new();

        while result == None {
            mypyself.queueimpl.unpark();
            backoff.snooze();

            // Check for exhaustion (or shutdown)...
            if mypyself.queueimpl.is_finished() {
                return Ok(None);
            }

            result = mypyself.queueimpl.queue.pop();
        }

        let result = result.unwrap();
        let gil = Python::acquire_gil();
        let py = gil.python();
        let pyout = PyDict::new(py);
        pyout
            .set_item("kmers", result.0)
            .expect("Error with Python");
        pyout
            .set_item("triple", result.1)
            .expect("Error with Python");

        // One last unparking...
        mypyself.queueimpl.unpark();

        Ok(Some(pyout.to_object(py)))
    }
}

/// Queue Impl
struct QueueImpl<Q> {
    // iter: Box<dyn Iterator<Item = I> + Send>,
    handles: Vec<JoinHandle<()>>,
    shutdown: Arc<AtomicBool>,
    exhausted: Arc<AtomicBool>,
    pub queue: Arc<ArrayQueue<Q>>,
}

impl<Q> QueueImpl<Q> {
    // fn new(iter: Box<dyn Iterator<Item = I> + Send>) -> Self {
    fn new<F>(queue_size: usize, threads: usize, func: F) -> Self
    where
        F: Fn(Arc<AtomicBool>, Arc<AtomicBool>, Arc<ArrayQueue<Q>>) + Sync + 'static + Send,
        Q: Send + 'static + Sync,
    {
        let shutdown = Arc::new(AtomicBool::new(false));
        let exhausted = Arc::new(AtomicBool::new(false));
        let queue = Arc::new(ArrayQueue::new(queue_size));

        let mut handles = Vec::new();

        for _i in 0..threads {

            let shutdown_c = Arc::clone(&shutdown);
            let exhausted_c = Arc::clone(&exhausted);
            let queue_c = Arc::clone(&queue);

            let handle = thread::spawn(move || {
                func(
                    Arc::clone(&shutdown_c),
                    Arc::clone(&exhausted_c),
                    Arc::clone(&queue_c),
                );
                exhausted_c.store(true, Ordering::SeqCst);
                shutdown_c.store(true, Ordering::SeqCst);
            });
            
            handles.push(handle);

        }

        
        QueueImpl {
            shutdown,
            exhausted,
            queue,
            handles,
        }
    }

    #[inline]
    fn is_finished(&self) -> bool {
        if self.queue.len() == 0
            && (self.exhausted.load(Ordering::Relaxed) || self.shutdown.load(Ordering::Relaxed))
        {
            return true;
        }
        return false;
    }

    #[inline]
    fn unpark(&self) {
        for i in &self.handles {
            i.thread().unpark();
        }
    }
}

// SequenceOrderKmersGenerator
type SequenceKmers = (WindowKmers, WindowKmers);
type SequenceOrder = bool;

type SequenceKmersSubmission = (SequenceKmers, SequenceOrder);

#[pyclass]
struct SequenceOrderKmersGenerator {
    wrapper: JoinHandle<()>,
    shutdown: Arc<AtomicBool>,
    exhausted: Arc<AtomicBool>,
    queue: Arc<ArrayQueue<SequenceKmersSubmission>>,
}

#[pymethods]
impl SequenceOrderKmersGenerator {
    #[new]
    fn new(k: usize, filename: String, window_size: usize, rand: bool, queue_size: usize) -> Self {
        let queue = Arc::new(ArrayQueue::new(queue_size));

        // Create KmerWindowGenerator
        let kmer_window_generator =
            KmerWindowGenerator::new(&filename, k, window_size * 2, 0, false, rand);

        let shutdown = Arc::new(AtomicBool::new(false));
        let exhausted = Arc::new(AtomicBool::new(false));

        let iter = Box::new(kmer_window_generator);

        let shutdown_c = Arc::clone(&shutdown);
        let exhausted_c = Arc::clone(&exhausted);

        let queue_handle = Arc::clone(&queue);
        let handle = thread::spawn(move || {
            let mut iter: Box<dyn Iterator<Item = KmerWindow> + Send> = iter;
            let mut offset = 0;
            let mut rc = false;

            let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);

            loop {
                let mut finished = false;

                while !finished {
                    let item = match iter.next() {
                        Some(x) => x,
                        None => {
                            offset += 1;

                            if (k == offset) && rc {
                                exhausted_c.store(true, Ordering::SeqCst);
                                shutdown_c.store(true, Ordering::SeqCst);
                                return;
                            } else {
                                if k == offset {
                                    rc = true;
                                    offset = 0;
                                }

                                let kmer_window_generator = KmerWindowGenerator::new(
                                    &filename,
                                    k,
                                    window_size,
                                    offset,
                                    rc,
                                    rand,
                                );

                                iter = Box::new(kmer_window_generator);
                                continue;
                            }
                        }
                    };

                    let KmerWindow { kmers, id, rc } = item;

                    if kmers.len() < window_size * 2 {
                        continue;
                    }

                    let kmers: Vec<Vec<bool>> = kmers
                        .iter()
                        .map(|x| convert_string_to_array(k, x))
                        .collect();
                    let first = kmers[0..window_size].to_vec();
                    let second = kmers[window_size..window_size * 2].to_vec();

                    let sequence_kmers;
                    let sequence_order;

                    if rng.gen::<bool>() {
                        sequence_kmers = (first, second);
                        sequence_order = false; // False is in order
                    } else {
                        sequence_kmers = (second, first);
                        sequence_order = true; // Out of order
                    }

                    let mut batch = (sequence_kmers, sequence_order);

                    while let Err(x) = queue_handle.push(batch) {
                        if shutdown_c.load(Ordering::Relaxed) {
                            // Mark as exhausted too then.
                            exhausted_c.store(true, Ordering::SeqCst);
                            return; // We are done, something triggered a
                                    // shutdown...
                        }
                        batch = x;
                        park(); // Queue is full, park the thread...
                    }
                }
            }
        });

        SequenceOrderKmersGenerator {
            wrapper: handle,
            shutdown,
            exhausted,
            queue,
        }
    }

    fn len(mypyself: PyRef<Self>) -> usize {
        mypyself.queue.len()
    }
}

#[pyproto]
impl PyIterProtocol for SequenceOrderKmersGenerator {
    fn __iter__(mypyself: PyRefMut<Self>) -> PyResult<PyObject> {
        let gil = Python::acquire_gil();
        let py = gil.python();
        Ok(mypyself.into_py(py))
    }

    fn __next__(mut mypyself: PyRefMut<Self>) -> PyResult<Option<PyObject>> {
        if mypyself.queue.len() == 0
            && (mypyself.exhausted.load(Ordering::Relaxed)
                || mypyself.shutdown.load(Ordering::Relaxed))
        {
            return Ok(None);
        }

        // Unpark the thread...
        mypyself.wrapper.thread().unpark();

        let mut result = mypyself.queue.pop();
        let backoff = Backoff::new();

        while result == None {
            mypyself.wrapper.thread().unpark();
            backoff.snooze();

            // Check for exhaustion (or shutdown)...
            if mypyself.exhausted.load(Ordering::Relaxed)
                || mypyself.shutdown.load(Ordering::Relaxed)
            {
                return Ok(None);
            }

            result = mypyself.queue.pop();
        }

        let result = result.unwrap();
        let gil = Python::acquire_gil();
        let py = gil.python();
        let pyout = PyDict::new(py);
        pyout
            .set_item("kmers", result.0)
            .expect("Error with Python");
        pyout
            .set_item("order", result.1)
            .expect("Error with Python");

        // One last unparking...
        mypyself.wrapper.thread().unpark();

        Ok(Some(pyout.to_object(py)))
    }
}

// ** Discriminator Masked Generator Wrapper
#[pyclass]
struct DiscriminatorMaskedGeneratorWrapper {
    iter: Box<dyn Iterator<Item = DiscriminatorMasked> + Send>, // + Send>,
    batch_size: usize,
    k: usize,
    offset: usize,
    window_size: usize,
    filename: String,
    replacement_pct: f32,
    rc: bool,
    rand: bool,
}

#[pyproto]
impl PyIterProtocol for DiscriminatorMaskedGeneratorWrapper {
    fn __iter__(mypyself: PyRefMut<Self>) -> PyResult<PyObject> {
        let gil = Python::acquire_gil();
        let py = gil.python();
        Ok(mypyself.into_py(py))
    }

    fn __next__(mut mypyself: PyRefMut<Self>) -> PyResult<Option<PyObject>> {
        // Generate Batch
        let mut batch_kmers: Vec<Vec<Vec<bool>>> = Vec::with_capacity(mypyself.batch_size);
        let mut batch_id: Vec<String> = Vec::with_capacity(mypyself.batch_size);
        let mut batch_truth: Vec<Vec<bool>> = Vec::with_capacity(mypyself.batch_size);

        while batch_kmers.len() < mypyself.batch_size {
            let item = match mypyself.iter.next() {
                Some(x) => x,
                None => {
                    mypyself.offset += 1;

                    if (mypyself.k == mypyself.offset) && mypyself.rc {
                        return Ok(None);
                    } else {
                        if mypyself.k == mypyself.offset {
                            mypyself.rc = true;
                            mypyself.offset = 0;
                        }

                        let kmer_window_generator = KmerWindowGenerator::new(
                            &mypyself.filename,
                            mypyself.k,
                            mypyself.window_size,
                            mypyself.offset,
                            mypyself.rc,
                            mypyself.rand,
                        );

                        let discriminator_masked_generator = DiscriminatorMaskedGenerator::new(
                            mypyself.replacement_pct,
                            mypyself.k,
                            kmer_window_generator,
                        );

                        mypyself.iter = Box::new(discriminator_masked_generator);
                        continue;
                    }
                }
            };

            let DiscriminatorMasked { kmers, id, truth } = item;
            let kmers = kmers
                .iter()
                .map(|x| convert_string_to_array(mypyself.k, x))
                .collect();
            batch_kmers.push(kmers);
            batch_id.push(id);
            batch_truth.push(truth);
        }

        let gil = Python::acquire_gil();
        let py = gil.python();
        let pyout = PyDict::new(py);
        pyout
            .set_item("kmers", batch_kmers)
            .expect("Error with Python");
        pyout.set_item("id", batch_id).expect("Error with Python");
        pyout
            .set_item("truth", batch_truth)
            .expect("Error with Python");
        Ok(Some(pyout.to_object(py)))
    }
}

#[pymethods]
impl DiscriminatorMaskedGeneratorWrapper {
    #[new]
    fn new(
        k: usize,
        filename: String,
        window_size: usize,
        batch_size: usize,
        replacement_pct: f32,
        rand: bool,
    ) -> Self {
        // Create KmerWindowGenerator
        let kmer_window_generator =
            KmerWindowGenerator::new(&filename, k, window_size, 0, false, rand);

        let discriminator_masked_generator =
            DiscriminatorMaskedGenerator::new(replacement_pct, k, kmer_window_generator);

        DiscriminatorMaskedGeneratorWrapper {
            iter: Box::new(discriminator_masked_generator),
            batch_size,
            k,
            offset: 0,
            filename,
            window_size,
            replacement_pct,
            rc: false,
            rand,
        }
    }
}

// Non-batch discriminator masked generator
#[pyclass]
struct DiscriminatorMaskedGeneratorWrapperNB {
    iter: Box<dyn Iterator<Item = DiscriminatorMasked> + Send>, // + Send>,
    k: usize,
    offset: usize,
    window_size: usize,
    filename: String,
    replacement_pct: f32,
    rc: bool,
}

#[pyproto]
impl PyIterProtocol for DiscriminatorMaskedGeneratorWrapperNB {
    fn __iter__(mypyself: PyRefMut<Self>) -> PyResult<PyObject> {
        let gil = Python::acquire_gil();
        let py = gil.python();
        Ok(mypyself.into_py(py))
    }

    fn __next__(mut mypyself: PyRefMut<Self>) -> PyResult<Option<PyObject>> {
        let mut finished = false;
        let mut item = None;

        while !finished {
            item = match mypyself.iter.next() {
                Some(x) => {
                    finished = true;
                    Some(x)
                }
                None => {
                    mypyself.offset += 1;
                    if (mypyself.k == mypyself.offset) && mypyself.rc {
                        println!("Finished, at the correct step...");
                        return Ok(None);
                    } else {
                        if mypyself.k == mypyself.offset {
                            mypyself.rc = true;
                            mypyself.offset = 0;
                        }

                        // println!("New Offset: {} {}", mypyself.offset, mypyself.rc);

                        let kmer_window_generator = KmerWindowGenerator::new(
                            &mypyself.filename,
                            mypyself.k,
                            mypyself.window_size,
                            mypyself.offset,
                            mypyself.rc,
                            true,
                        );

                        let discriminator_masked_generator = DiscriminatorMaskedGenerator::new(
                            mypyself.replacement_pct,
                            mypyself.k,
                            kmer_window_generator,
                        );

                        mypyself.iter = Box::new(discriminator_masked_generator);
                        continue;
                    }
                }
            };
        }

        match item {
            Some(x) => {
                let DiscriminatorMasked {
                    kmers,
                    truth,
                    id: _,
                } = x;
                let kmers: Vec<Vec<bool>> = kmers
                    .iter()
                    .map(|x| convert_string_to_array(mypyself.k, x))
                    .collect();

                let gil = Python::acquire_gil();
                let py = gil.python();
                let pyout = PyDict::new(py);
                // let pyout = PyTuple::new(py, [kmers, truth]);
                pyout.set_item("kmers", kmers).expect("Py Error");
                pyout.set_item("truths", truth).expect("Py Error");
                Ok(Some(pyout.to_object(py)))
            }
            None => Ok(None),
        }
    }
}

#[pymethods]
impl DiscriminatorMaskedGeneratorWrapperNB {
    #[new]
    fn new(k: usize, filename: String, window_size: usize, replacement_pct: f32) -> Self {
        // Create KmerWindowGenerator
        let kmer_window_generator =
            KmerWindowGenerator::new(&filename, k, window_size, 0, false, true);

        let discriminator_masked_generator =
            DiscriminatorMaskedGenerator::new(replacement_pct, k, kmer_window_generator);

        DiscriminatorMaskedGeneratorWrapperNB {
            iter: Box::new(discriminator_masked_generator),
            k,
            offset: 0,
            filename,
            window_size,
            replacement_pct,
            rc: false,
        }
    }
}

#[pyfunction]
fn convert_fasta_to_sfasta(input: String, output: String) {
    sfasta::convert_fasta_file(&input, &output);
}

#[pyfunction]
fn index_sfasta(input: String) -> String {
    sfasta::index(&input)
}

#[pyfunction]
fn get_headers_from_sfasta(input: String) -> Vec<String> {
    sfasta::get_headers_from_sfasta(input)
}

#[pyfunction]
fn test_sfasta(input: String) {
    sfasta::test_sfasta(input);
}

// ** Fasta Kmer Generator
#[pyclass]
struct FastaKmersGenerator {
    iter: Box<dyn Iterator<Item = KmerCoordsWindow> + Send>,
    k: usize,
    offset: usize,
    window_size: usize,
    filename: String,
    rc: bool,
    sliding: bool,
    start_rc: bool,
}

#[pyproto]
impl PyIterProtocol for FastaKmersGenerator {
    fn __iter__(mypyself: PyRefMut<Self>) -> PyResult<PyObject> {
        let gil = Python::acquire_gil();
        let py = gil.python();
        Ok(mypyself.into_py(py))
    }

    fn __next__(mut mypyself: PyRefMut<Self>) -> PyResult<Option<PyObject>> {
        let mut finished = false;
        let mut item = None;

        while !finished {
            item = match mypyself.iter.next() {
                Some(x) => {
                    finished = true;
                    Some(x)
                }
                None => {
                    mypyself.offset += 1;
                    if (mypyself.k == mypyself.offset && mypyself.rc) || !mypyself.sliding {
                        return Ok(None);
                    } else {
                        if mypyself.k == mypyself.offset {
                            mypyself.rc = true;
                            mypyself.offset = 0;
                        }

                        let iter = KmerCoordsWindowIter::new(
                            &mypyself.filename,
                            mypyself.k,
                            mypyself.window_size,
                            mypyself.offset,
                            mypyself.rc,
                            false,
                        );

                        mypyself.iter = Box::new(iter);
                        continue;
                    }
                }
            };
        }

        match item {
            Some(x) => {
                let KmerCoordsWindow {
                    kmers,
                    coords,
                    id,
                    rc,
                } = x;
                let kmers: Vec<Vec<bool>> = kmers
                    .iter()
                    .map(|x| convert_string_to_array(mypyself.k, x))
                    .collect();

                let gil = Python::acquire_gil();
                let py = gil.python();
                let pyout = PyDict::new(py);
                // let pyout = PyTuple::new(py, [kmers, truth]);
                pyout.set_item("kmers", kmers).expect("Py Error");
                pyout.set_item("coords", coords).expect("Py Error");
                pyout.set_item("ids", id).expect("Py Error");
                pyout.set_item("rc", rc).expect("Py Error");
                Ok(Some(pyout.to_object(py)))
            }
            None => Ok(None),
        }
    }
}

#[pymethods]
impl FastaKmersGenerator {
    /// Create a new FastaKmersGenerator
    /// Arguments:
    ///     k: kmer size, integer (usize)
    ///     filename: String, location of .fasta, .sfata file
    ///     window_size: How many kmers to produce per iteration
    ///     sliding: Sliding window or just a once-over?
    ///     start_rc: Start on the reverse strand? WARNING: Only use this when
    /// not doing sliding windows...
    /// TODO: Need to give this a chance to return less than window_size kmer's
    #[new]
    fn new(k: usize, filename: String, window_size: usize, sliding: bool, start_rc: bool) -> Self {
        // Create KmerWindowGenerator
        let iter = KmerCoordsWindowIter::new(&filename.clone(), k, window_size, 0, start_rc, false);

        FastaKmersGenerator {
            iter: Box::new(iter),
            k,
            offset: 0,
            filename,
            window_size,
            rc: false,
            sliding,
            start_rc,
        }
    }
}

/// Provides functions for python dealing with Kmers from fasta and sfasta
/// files...
#[pymodule]
fn pyracular(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<MaskedKmersGenerator>()?;
    m.add_class::<DiscriminatorMaskedGeneratorWrapper>()?;
    m.add_class::<DiscriminatorMaskedGeneratorWrapperNB>()?;
    m.add_class::<Gff3KmerGenerator>()?;
    m.add_class::<FastaKmersGenerator>()?;
    m.add_class::<SequenceOrderKmersGenerator>()?;
    m.add_class::<MatchedKmersGenerator>()?;
    m.add_class::<TripleLossKmersGenerator>()?;
    m.add_wrapped(wrap_pyfunction!(convert_fasta_to_sfasta))?;
    m.add_wrapped(wrap_pyfunction!(get_headers_from_sfasta))?;
    m.add_wrapped(wrap_pyfunction!(test_sfasta))?;
    m.add_wrapped(wrap_pyfunction!(index_sfasta))?;

    Ok(())
}
