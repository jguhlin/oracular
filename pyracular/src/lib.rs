use crossbeam::queue::ArrayQueue;
use crossbeam::utils::Backoff;
use std::sync::atomic::Ordering;
use std::thread;
use std::thread::{park, JoinHandle};

use rand::prelude::*;
use std::sync::atomic::AtomicBool;
use std::sync::Arc;

use mimalloc::MiMalloc;
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use rand_xoshiro::rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

// NOTE: New naming convention
// Rust-y stuff is "iter" Python is "Generator"

use liboracular::io;
use liboracular::kmers::KmerWindow;
use liboracular::kmers::{
    rc_kmerwindow, replace_blank, replace_random, DiscriminatorMasked,
    DiscriminatorMaskedGenerator, Gff3Kmers, Gff3KmersIter, KmerCoordsWindow, KmerCoordsWindowIter,
    KmerWindowGenerator,
};
use libsfasta::prelude::*;

use pyo3::class::iter::IterNextOutput;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use pyo3::wrap_pyfunction;

pub mod taxonomy;
pub use taxonomy::*;

#[inline]
fn convert_string_to_array(k: usize, s: &[u8]) -> Vec<bool> {
    let mut out: Vec<bool> = vec![false; k * 5];

    for (x, c) in s.iter().enumerate() {
        match c {
            65 => out[5 * x] = true,     // A
            84 => out[5 * x + 1] = true, // T
            67 => out[5 * x + 3] = true, // C
            71 => out[5 * x + 4] = true, // G
            78 => out[5 * x + 2] = true, // N
            _ => out[5 * x + 2] = true,  // N for everything else...
        };
    }

    out
}

#[pyfunction]
// The end is chopped off.
// TODO: Handle if longer than k + ws
fn convert_sequence_to_array(k: usize, ws: usize, s: &str) -> Vec<bool> {
    let mut out: Vec<bool> = Vec::with_capacity(s.len() * 5);

    let kmers = f32::floor(s.len() as f32 / k as f32) as usize;
    for i in 0..kmers {
        let mut kmer = convert_string_to_array(k, &s[i * k..(i + 1) * k].as_bytes());
        out.append(&mut kmer);
    }

    out
}

#[pyfunction]
// Given k as kmer length
// and window_size as the window size
// Generate a mask for each kmer (1 to pay attention to it, 0 to ignore)
// Return is the new sequence (with padding) and the mask
fn pad_and_mask(k: usize, window_size: usize, mut seq: Vec<bool>) -> (Vec<bool>, Vec<bool>) {
    let mut mask: Vec<bool> = vec![false; window_size];

    for i in 0..(seq.len() / (5 * k)) as usize {
        mask[i] = true;
    }

    if seq.len() < window_size * k * 5 {
        let mut pad = vec![false; window_size * k * 5 - seq.len()];
        seq.append(&mut pad);
    }

    (seq, mask)
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

#[pymethods]
impl Gff3KmerGenerator {
    fn __iter__(mypyself: PyRef<Self>) -> PyRef<Self> {
        mypyself
    }

    fn __next__(mut mypyself: PyRefMut<Self>) -> IterNextOutput<PyObject, &'static str> {
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
                        return IterNextOutput::Return("Finished");
                    } else {
                        if mypyself.k == mypyself.offset {
                            mypyself.rc = true;
                            mypyself.offset = 0;
                        }

                        let kmercoords_window_iter = KmerCoordsWindowIter::new(
                            mypyself.filename.clone(),
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

                Python::with_gil(|py| -> IterNextOutput<pyo3::Py<pyo3::PyAny>, &str> {
                    let result = (kmers, classifications, coords, id, rc);
                    IterNextOutput::Yield(result.to_object(py))
                })
            }
            None => IterNextOutput::Return("Finished"),
        }
    }

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
            KmerCoordsWindowIter::new(filename.clone(), k, window_size, 0, false, rand);

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
        seed: usize,
    ) -> Self {
        let queueimpl = QueueImpl::new(
            queue_size,
            24,
            seed as u64,
            move |shutdown, _exhausted, queue, seed, thread_number| {
                let mut offset = 0;
                let mut rc = false;

                // TODO: seed and thread number?

                loop {
                    // Create KmerWindowGenerator
                    let kmer_window_generator = KmerWindowGenerator::new(
                        filename.clone(),
                        k,
                        window_size,
                        offset,
                        rc,
                        rand,
                    );

                    let discriminator_masked_generator = DiscriminatorMaskedGenerator::new(
                        replacement_pct,
                        k,
                        kmer_window_generator,
                    );

                    let iter = Box::new(discriminator_masked_generator);

                    for x in iter {
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
            },
        );

        MaskedKmersGenerator { queueimpl }
    }

    fn len(mypyself: PyRef<Self>) -> usize {
        mypyself.queueimpl.queue.len()
    }

    fn __iter__(mypyself: PyRefMut<Self>) -> PyResult<PyObject> {
        Python::with_gil(|py| -> PyResult<PyObject> { Ok(mypyself.into_py(py)) })
    }

    fn __next__(mypyself: PyRef<Self>) -> PyResult<Option<PyObject>> {
        if mypyself.queueimpl.is_finished() {
            return Ok(None);
        }

        // Unpark the thread...
        mypyself.queueimpl.unpark();

        let mut result = mypyself.queueimpl.queue.pop();
        let backoff = Backoff::new();

        while result.is_none() {
            mypyself.queueimpl.unpark();
            backoff.snooze();

            // Check for exhaustion (or shutdown)...
            if mypyself.queueimpl.is_finished() {
                return Ok(None);
            }

            result = mypyself.queueimpl.queue.pop();
        }

        let result = result.unwrap();
        Python::with_gil(|py| -> PyResult<Option<PyObject>> {
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
        })
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
    fn new(k: usize, filename: String, window_size: usize, queue_size: usize, seed: usize) -> Self {
        let queueimpl = QueueImpl::new(
            queue_size,
            16,
            seed as u64,
            move |shutdown, _exhausted, queue, seed, thread_number| {
                let mut offset = 0;
                let mut rc = false;
                let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);

                for i in 0..thread_number {
                    rng.jump();
                }

                // TODO: Make even smarter -- Load up 1k windows and pick from there matching
                // and non-matching ones, including some RC ones as well...
                loop {
                    // Create KmerWindowGenerator
                    let mut iter1 = KmerWindowGenerator::new(
                        filename.clone(),
                        k,
                        window_size,
                        offset,
                        rc,
                        true,
                    );

                    let mut iter2 = KmerWindowGenerator::new(
                        filename.clone(),
                        k,
                        window_size,
                        offset,
                        rc,
                        true,
                    );

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
            },
        );

        MatchedKmersGenerator { queueimpl }
    }

    fn len(mypyself: PyRef<Self>) -> usize {
        mypyself.queueimpl.queue.len()
    }

    fn __iter__(mypyself: PyRefMut<Self>) -> PyResult<PyObject> {
        Python::with_gil(|py| -> PyResult<PyObject> { Ok(mypyself.into_py(py)) })
    }

    fn __next__(mypyself: PyRef<Self>) -> PyResult<Option<PyObject>> {
        if mypyself.queueimpl.is_finished() {
            return Ok(None);
        }

        // Unpark the thread...
        mypyself.queueimpl.unpark();

        let mut result = mypyself.queueimpl.queue.pop();
        let backoff = Backoff::new();

        while result.is_none() {
            mypyself.queueimpl.unpark();
            backoff.snooze();

            // Check for exhaustion (or shutdown)...
            if mypyself.queueimpl.is_finished() {
                return Ok(None);
            }

            result = mypyself.queueimpl.queue.pop();
        }

        let result = result.unwrap();
        Python::with_gil(|py| -> PyResult<_> {
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
        })
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
    queueimpl: QueueImpl<TripleLossReturn>,
}

type ReverseComplement = bool;
type Truths = Vec<bool>;

type TripleLossSubmission = (MatchedKmers, (Truths, Truths, Matches, ReverseComplement));

// TODO: Multiple threads is good, make it actually divide the work instead of
// just repeating the work...

#[pymethods]
impl TripleLossKmersGenerator {
    #[new]
    fn new(
        k: usize,
        filename: String,
        replacement_pct: f32,
        window_size: usize,
        queue_size: usize,
        threads: usize,
        seed: usize,
    ) -> Self {
        let queueimpl = QueueImpl::new(
            queue_size,
            threads,
            seed as u64,
            move |shutdown, exhausted, queue, seed, thread_number| {
                let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
                for _i in 0..thread_number {
                    rng.jump();
                }

                // let mut sfasta = SfastaParser::open(filename.clone()).expect("Unable to open
                // file"); TODO: Make sfasta shareable across threads
                let mut buf = std::io::BufReader::new(
                    std::fs::File::open(filename.clone()).expect("Unable to open file"),
                );
                let mut sfasta = SfastaParser::open_from_buffer(&mut buf, false).unwrap();
                sfasta.get_seqlocs().expect("Unable to get seqlocs");

                let mut valid_indices = Vec::new();
                for i in 0..sfasta.seqlocs.as_ref().unwrap().total_seqlocs {
                    let seqloc = sfasta.get_seqloc(i).expect("Unable to get seqloc");

                    if let Some(seqloc) = seqloc {
                        if seqloc.sequence.is_some()
                            && sfasta.seqloc_len(&seqloc) >= k * window_size + k
                        {
                            valid_indices.push(i);
                        }
                    }
                }

                log::debug!("Generating kmers for thread {}", thread_number);

                // Disable masking
                sfasta.masking = None;

                let minimum_seqlength = k * window_size + k;

                let total_seqlocs = sfasta.seqlocs.as_ref().unwrap().total_seqlocs;

                // TODO: Make even smarter -- Load up 1k windows and pick from there matching
                // and non-matching ones, including some RC ones as well...

                let choice_dist = rand::distributions::WeightedIndex::new([1, 1, 2]).unwrap();
                let offset_dist = rand::distributions::Uniform::from(0..k);
                let seqloc_dist = rand::distributions::Uniform::from(0..valid_indices.len());

                loop {
                    // let offset = rng.gen_range(0..k);
                    let offset = offset_dist.sample(&mut rng);
                    let rc: bool = rng.gen();
                    // let seqloc_i = rng.gen_range(0..total_seqlocs);
                    let seqloc_i = valid_indices[seqloc_dist.sample(&mut rng)];
                    let seqloc = sfasta.get_seqloc(seqloc_i).unwrap().unwrap();

                    let sequence = match sfasta.get_sequence_only_by_seqloc(&seqloc, true) {
                        Ok(Some(x)) => x,
                        Ok(None) => panic!("Unable to get sequence"),
                        Err(_) => continue,
                    };

                    if sequence.len() < minimum_seqlength {
                        log::debug!("Skipping sequence due to length");
                        continue;
                    }

                    let mut iter1 =
                    // KmerWindowGenerator::new(&filename, k, window_size, offset, rc, true);
                    KmerWindowGenerator::from_sequence(sequence, k, window_size, offset, rc).unwrap();

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
                        // Non-matched sequence -- Kmer window 1 and 2 from completely different
                        // seqs...
                        // let choice: u8 = rng.gen_range(0..3); // Give us a number between 0 and 2
                        let choice: u8 = choice_dist.sample(&mut rng) as u8;

                        // Always need a starting window...
                        item1 = match iter1.next() {
                            Some(x) => x,
                            None => {
                                break 'inner;
                            }
                        };

                        while item1.kmers.len() < window_size {
                            log::debug!("{} {}", item1.kmers.len(), window_size);
                            item1 = match iter1.next() {
                                Some(x) => {
                                    if x.kmers.is_empty() {
                                        break 'inner;
                                    }
                                    x
                                }
                                None => {
                                    break 'inner;
                                }
                            };
                        }

                        let start = std::time::Instant::now();

                        // Matched Sequence
                        if choice == 0 {
                            matched = true;
                            reversecomplement = false;

                            item2 = match get_random_sequence_from_seqloc(
                                &mut sfasta,
                                k,
                                window_size,
                                &seqloc,
                                &mut rng,
                                true,
                            ) {
                                Some(x) => x,
                                None => continue,
                            };

                            log::debug!("Choice 0: {:#?}", start.elapsed());

                        // RC
                        } else if choice == 1 {
                            matched = true;
                            reversecomplement = true;

                            item2 = item1.clone();
                            item2 = rc_kmerwindow(item2);

                            log::debug!("Choice 1: {:#?}", start.elapsed());
                        // Not matching sequence...
                        } else {
                            matched = false;
                            reversecomplement = false;

                            // let mut seqloc_j = rng.gen_range(0..total_seqlocs);
                            let mut seqloc_j = valid_indices[seqloc_dist.sample(&mut rng)];

                            // Unlikely, but to be sure...
                            while seqloc_i == seqloc_j {
                                //seqloc_j = rng.gen_range(0..total_seqlocs);
                                seqloc_j = valid_indices[seqloc_dist.sample(&mut rng)];
                            }

                            let rand_seqloc = sfasta.get_seqloc(seqloc_j).unwrap().unwrap();

                            item2 = match get_random_sequence_from_seqloc(
                                &mut sfasta,
                                k,
                                window_size,
                                &rand_seqloc,
                                &mut rng,
                                false,
                            ) {
                                Some(x) => x,
                                None => continue,
                            };

                            log::debug!("Choice 2: {:#?}", start.elapsed());
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

                        // kmers3
                        let kmers1_truth = kmers1
                            .clone()
                            .iter()
                            .map(|x| convert_string_to_array(k, x))
                            .collect();

                        let kmers2_truth = kmers2
                            .clone()
                            .iter()
                            .map(|x| convert_string_to_array(k, x))
                            .collect();

                        let truth1 = replace_blank(k, replacement_pct, &mut kmers1, &mut rng);
                        let truth2 = replace_blank(k, replacement_pct, &mut kmers2, &mut rng);

                        let kmers1 = kmers1
                            .iter()
                            .map(|x| convert_string_to_array(k, x))
                            .collect();

                        let kmers2 = kmers2
                            .iter()
                            .map(|x| convert_string_to_array(k, x))
                            .collect();

                        let mut batch = TripleLossReturn {
                            kmers1,
                            kmers2,
                            kmers3: kmers1_truth,
                            kmers4: kmers2_truth,
                            truth1,
                            truth2,
                            matched,
                            reversecomplement,
                        };

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
                }
            },
        );

        TripleLossKmersGenerator { queueimpl }
    }

    fn len(mypyself: PyRef<Self>) -> usize {
        mypyself.queueimpl.queue.len()
    }
    fn __iter__(mypyself: PyRef<Self>) -> PyRef<Self> {
        mypyself
    }

    fn __next__(mut mypyself: PyRefMut<Self>) -> IterNextOutput<TripleLossReturn, &'static str> {
        let mut result = mypyself.queueimpl.queue.pop();
        let backoff = Backoff::new();

        // TODO: python allow_threads
        while result.is_none() {
            // Unpark the threads...
            mypyself.queueimpl.unpark();
            backoff.snooze();

            // Check for exhaustion (or shutdown)...
            if mypyself.queueimpl.is_finished() {
                mypyself.queueimpl.shutdown();
                return IterNextOutput::Return("Finished");
            }

            result = mypyself.queueimpl.queue.pop();
        }

        let result = result.unwrap();

        IterNextOutput::Yield(result)
    }
}

#[pyclass(frozen)]
pub struct TripleLossReturn {
    #[pyo3(get)]
    pub kmers1: Vec<Vec<bool>>,
    #[pyo3(get)]
    pub kmers2: Vec<Vec<bool>>,
    #[pyo3(get)]
    pub kmers3: Vec<Vec<bool>>,
    #[pyo3(get)]
    pub kmers4: Vec<Vec<bool>>,
    #[pyo3(get)]
    pub truth1: Vec<bool>,
    #[pyo3(get)]
    pub truth2: Vec<bool>,
    #[pyo3(get)]
    pub matched: bool,
    #[pyo3(get)]
    pub reversecomplement: bool,
}

impl Drop for TripleLossKmersGenerator {
    fn drop(&mut self) {
        self.queueimpl.shutdown();
    }
}

/// Support functions for triple loss generator
fn is_all_ns(seq: &[u8]) -> bool {
    bytecount::count(seq, b'N') == seq.len()
}

/// Support functions for triple loss generator
fn get_random_sequence_from_id<R: Rng + ?Sized>(
    sfasta: &mut Sfasta,
    k: usize,
    window_size: usize,
    id: &str,
    rng: &mut R,
) -> Option<KmerWindow> {
    let needed_length = (k * window_size) + k;

    let mut seq;

    seq = sfasta.get_sequence_by_id(id).unwrap().unwrap();
    if seq.len() < needed_length {
        return None;
    }

    let seqlen = seq.len().saturating_sub(needed_length);

    if seqlen == 0 {
        return None;
    }

    let mut start = rng.gen_range(0..seqlen);
    let mut end = start + needed_length;
    assert!(end < seq.len());

    while is_all_ns(&seq.sequence.as_ref().unwrap()[start..end]) {
        start = rng.gen_range(0..seqlen);
        end = start + needed_length;
    }

    seq.sequence = Some(seq.sequence.unwrap()[start..end].to_vec());

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
fn get_random_sequence_from_seqloc<R: Rng + ?Sized>(
    sfasta: &mut Sfasta,
    k: usize,
    window_size: usize,
    seqloc: &libsfasta::datatypes::SeqLoc,
    rng: &mut R,
    caching: bool,
) -> Option<KmerWindow> {
    let needed_length = (k * window_size) + k;

    let mut seq;

    let seqlen = sfasta.seqloc_len(seqloc) as u32;

    let needed_length = needed_length as u32;

    let seqlen = seqlen.saturating_sub(needed_length);

    if seqlen == 0 {
        return None;
    }

    let mut start = rng.gen_range(0..seqlen);
    let mut end = start + needed_length;

    let locs = sfasta.seq_slice(seqloc, start..end);
    seq = sfasta
        .get_sequence_only_by_locs(&locs, caching)
        .unwrap()
        .unwrap();

    while is_all_ns(seq.sequence.as_ref().unwrap()) {
        start = rng.gen_range(0..seqlen);
        end = start + needed_length;
        let locs = sfasta.seq_slice(seqloc, start..end);
        seq = sfasta
            .get_sequence_only_by_locs(&locs, caching)
            .unwrap()
            .unwrap();
    }

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
    sfasta: &mut Sfasta,
    k: usize,
    window_size: usize,
    locs: &Vec<u64>,
    mut rng: &mut R,
) -> Option<KmerWindow> {
    let needed_length = ((k * window_size) + k) as u64;

    let mut loc = *locs.choose(&mut rng).unwrap();
    let mut orig_seq = sfasta
        .get_sequence_by_index(loc as usize)
        .expect("Unable to get seq")
        .expect("Unable to get seq");

    // TODO: Move length filter elsewhere
    while orig_seq.len() < needed_length as usize {
        loc = *locs.choose(&mut rng).unwrap();
        orig_seq = sfasta.get_sequence_by_index(loc as usize).unwrap().unwrap();
    }

    let seqlen = orig_seq.len().saturating_sub(needed_length as usize);

    if seqlen == 0 {
        return None;
    }

    let start = rng.gen_range(0..seqlen);
    let mut end = start + needed_length as usize;
    assert!(end < orig_seq.len());
    if orig_seq.len() >= end + 1000 {
        end += 1000;
    } else {
        end = orig_seq.len();
    }

    // TODO: Replace function (and above!)
    // So that we only extract the sequence we need (esp. for large seq's which may
    // be in multiple blocks)
    /*let sequence = sfasta
    .get_seq_slice(seq.id.clone(), start, end)
    .expect("Unable to get seq slice"); */

    let workseq = io::Sequence {
        sequence: Some(orig_seq.sequence.unwrap()[start..end].to_vec()),
        scores: None,
        header: None,
        id: orig_seq.id,
        offset: orig_seq.offset,
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

/// Queue Impl
struct QueueImpl<Q> {
    // iter: Box<dyn Iterator<Item = I> + Send>,
    handles: Vec<JoinHandle<()>>,
    shutdown: Arc<AtomicBool>,
    exhausted: Vec<Arc<AtomicBool>>,
    pub queue: Arc<ArrayQueue<Q>>,
}

impl<Q> QueueImpl<Q> {
    // fn new(iter: Box<dyn Iterator<Item = I> + Send>) -> Self {
    fn new<F>(queue_size: usize, threads: usize, seed: u64, func: F) -> Self
    where
        F: Fn(Arc<AtomicBool>, Arc<AtomicBool>, Arc<ArrayQueue<Q>>, u64, usize)
            + Sync
            + 'static
            + Send,
        Q: Send + 'static + Sync,
    {
        let shutdown = Arc::new(AtomicBool::new(false));
        let mut exhausted = Vec::new();
        let queue = Arc::new(ArrayQueue::new(queue_size));

        let mut handles = Vec::new();

        let fc = Arc::new(func);

        for i in 0..threads {
            let exhausted_handle = Arc::new(AtomicBool::new(false));
            let shutdown_c = Arc::clone(&shutdown);
            let exhausted_c = Arc::clone(&exhausted_handle);
            let queue_c = Arc::clone(&queue);

            exhausted.push(exhausted_handle);

            let f = Arc::clone(&fc);

            let handle = thread::spawn(move || {
                f(
                    Arc::clone(&shutdown_c),
                    Arc::clone(&exhausted_c),
                    Arc::clone(&queue_c),
                    seed,
                    i,
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
        let mut exhausted = true;
        // Are all exhausted?
        for handle in &self.exhausted {
            if !handle.load(Ordering::Relaxed) {
                exhausted = false;
            }
        }

        self.queue.len() == 0 && (exhausted || self.shutdown.load(Ordering::Relaxed))
    }

    #[inline]
    fn unpark(&self) {
        for i in &self.handles {
            i.thread().unpark();
        }
    }

    fn shutdown(&mut self) {
        self.shutdown.store(true, Ordering::SeqCst);
        self.unpark();
        for i in self.handles.drain(..) {
            i.join().expect("Unable to shut down");
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
            KmerWindowGenerator::new(filename.clone(), k, window_size * 2, 0, false, rand);

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
                                filename.clone(),
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

    fn __iter__(mypyself: PyRefMut<Self>) -> PyResult<PyObject> {
        Python::with_gil(|py| -> PyResult<_> { Ok(mypyself.into_py(py)) })
    }

    fn __next__(mypyself: PyRef<Self>) -> PyResult<Option<PyObject>> {
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
        Python::with_gil(|py| -> PyResult<_> {
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
        })
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

#[pymethods]
impl DiscriminatorMaskedGeneratorWrapper {
    fn __iter__(mypyself: PyRefMut<Self>) -> PyResult<PyObject> {
        Python::with_gil(|py| -> PyResult<_> { Ok(mypyself.into_py(py)) })
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
                            mypyself.filename.clone(),
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

        Python::with_gil(|py| -> PyResult<_> {
            let pyout = PyDict::new(py);
            pyout
                .set_item("kmers", batch_kmers)
                .expect("Error with Python");
            pyout.set_item("id", batch_id).expect("Error with Python");
            pyout
                .set_item("truth", batch_truth)
                .expect("Error with Python");
            Ok(Some(pyout.to_object(py)))
        })
    }

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
            KmerWindowGenerator::new(filename.clone(), k, window_size, 0, false, rand);

        let discriminator_masked_generator =
            DiscriminatorMaskedGenerator::new(replacement_pct, k, kmer_window_generator);

        DiscriminatorMaskedGeneratorWrapper {
            iter: Box::new(discriminator_masked_generator),
            batch_size,
            k,
            offset: 0,
            filename: filename.clone(),
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

#[pymethods]
impl DiscriminatorMaskedGeneratorWrapperNB {
    fn __iter__(mypyself: PyRefMut<Self>) -> PyResult<PyObject> {
        Python::with_gil(|py| -> PyResult<_> { Ok(mypyself.into_py(py)) })
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
                            mypyself.filename.clone(),
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

                Python::with_gil(|py| -> PyResult<_> {
                    let pyout = PyDict::new(py);
                    // let pyout = PyTuple::new(py, [kmers, truth]);
                    pyout.set_item("kmers", kmers).expect("Py Error");
                    pyout.set_item("truths", truth).expect("Py Error");
                    Ok(Some(pyout.to_object(py)))
                })
            }
            None => Ok(None),
        }
    }

    #[new]
    fn new(k: usize, filename: String, window_size: usize, replacement_pct: f32) -> Self {
        // Create KmerWindowGenerator
        let kmer_window_generator =
            KmerWindowGenerator::new(filename.clone(), k, window_size, 0, false, true);

        let discriminator_masked_generator =
            DiscriminatorMaskedGenerator::new(replacement_pct, k, kmer_window_generator);

        DiscriminatorMaskedGeneratorWrapperNB {
            iter: Box::new(discriminator_masked_generator),
            k,
            offset: 0,
            filename: filename.clone(),
            window_size,
            replacement_pct,
            rc: false,
        }
    }
}

#[pyfunction]
fn convert_fasta_to_sfasta(input: String, output: String) {
    convert_fasta_file(&input, &output);
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
    rand: bool,
}

#[pymethods]
impl FastaKmersGenerator {
    fn __iter__(mypyself: PyRef<Self>) -> PyRef<Self> {
        mypyself
    }

    fn __next__(mut mypyself: PyRefMut<Self>) -> IterNextOutput<PyObject, &'static str> {
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
                    println!(
                        "Starting a new one... {:#?} {:#?} {:#?} {:#?}",
                        mypyself.k, mypyself.offset, mypyself.rc, mypyself.sliding
                    );
                    if (mypyself.k == mypyself.offset && mypyself.rc) || !mypyself.sliding {
                        return IterNextOutput::Return("Finished");
                    } else {
                        if mypyself.k == mypyself.offset {
                            mypyself.rc = true;
                            mypyself.offset = 0;
                        }

                        let iter = KmerCoordsWindowIter::new(
                            mypyself.filename.clone(),
                            mypyself.k,
                            mypyself.window_size,
                            mypyself.offset,
                            mypyself.rc,
                            mypyself.rand,
                        );

                        mypyself.iter = Box::new(iter);
                        println!("Getting next...");
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

                Python::with_gil(|py| -> IterNextOutput<_, _> {
                    //                let pyout = PyDict::new(py);
                    // let pyout = PyTuple::new(py, [kmers, truth]);
                    //                pyout.set_item("kmers", kmers).expect("Py Error");
                    //                pyout.set_item("coords", coords).expect("Py Error");
                    //                pyout.set_item("ids", id).expect("Py Error");
                    //                pyout.set_item("rc", rc).expect("Py Error");
                    let result = (kmers, coords, id, rc);
                    return IterNextOutput::Yield(result.to_object(py));
                })
            }
            None => IterNextOutput::Return("Finished"),
        }
    }

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
    fn new(
        k: usize,
        filename: String,
        window_size: usize,
        sliding: bool,
        start_rc: bool,
        rand: bool,
    ) -> Self {
        // Create KmerWindowGenerator
        let iter = KmerCoordsWindowIter::new(filename.clone(), k, window_size, 0, start_rc, rand);

        FastaKmersGenerator {
            iter: Box::new(iter),
            k,
            offset: 0,
            filename,
            window_size,
            rc: false,
            sliding,
            start_rc,
            rand,
        }
    }
    /*
    fn get_next_seq(&mut self) {
        self.iter.get_next_seq();
    }*/
}

/// Provides functions for python dealing with Kmers from fasta and sfasta
/// files...
#[pymodule]
fn pyracular(_py: Python, m: &PyModule) -> PyResult<()> {
    env_logger::init();
    m.add_class::<TripleLossReturn>()?;
    m.add_class::<MaskedKmersGenerator>()?;
    m.add_class::<DiscriminatorMaskedGeneratorWrapper>()?;
    m.add_class::<DiscriminatorMaskedGeneratorWrapperNB>()?;
    m.add_class::<Gff3KmerGenerator>()?;
    m.add_class::<FastaKmersGenerator>()?;
    m.add_class::<SequenceOrderKmersGenerator>()?;
    m.add_class::<MatchedKmersGenerator>()?;
    m.add_class::<TripleLossKmersGenerator>()?;
    // m.add_wrapped(wrap_pyfunction!(convert_fasta_to_sfasta))?;
    m.add_wrapped(wrap_pyfunction!(pad_and_mask));
    m.add_wrapped(wrap_pyfunction!(convert_sequence_to_array))?;

    Ok(())
}
