extern crate crossbeam;
extern crate mimalloc;

use crossbeam::queue::ArrayQueue;
use crossbeam::utils::Backoff;
use mimalloc::MiMalloc;
use std::sync::atomic::Ordering;
use std::thread;
use std::thread::{park, JoinHandle};

use std::sync::atomic::AtomicBool;
use std::sync::Arc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

// NOTE: New naming convention
// Rust-y stuff is "iter" Python is "Generator"

use liboracular::kmers::{
    DiscriminatorMasked, DiscriminatorMaskedGenerator, Gff3Kmers, Gff3KmersIter, KmerCoordsWindow,
    KmerCoordsWindowIter, KmerWindowGenerator,
};
use liboracular::sfasta;

use pyo3::prelude::*;
// use pyo3::wrap_pyfunction;
use pyo3::types::PyDict;
use pyo3::PyIterProtocol;

// use pyo3::wrap_pyfunction;

use pyo3::wrap_pyfunction;

// use std::fs::File;
// use std::io::BufReader;

// TODO: Should be literals to concatenate...
#[inline(always)]
fn convert_string_to_array(k: usize, s: &[u8]) -> Vec<u8> {
    let mut out: Vec<u8> = vec![0; k * 5];

    for (x, c) in s.iter().enumerate() {
        match c {
            65 => out[5*x] = 1,   // A
            84 => out[5*x+1] = 1, // T
            67 => out[5*x+3] = 1, // C
            71 => out[5*x+4] = 1, // G
            78 => out[5*x+2] = 1, // N
            _  => out[5*x+2] = 1, // N for everything else...
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
                let kmers: Vec<Vec<u8>> = kmers
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
    wrapper: JoinHandle<()>,
    shutdown: Arc<AtomicBool>,
    exhausted: Arc<AtomicBool>,
    queue: Arc<ArrayQueue<WindowSubmission>>,
}

type KmerStr = String;
type Kmer = Vec<u8>;
type Id = String;
type Truth = u8;

type WindowKmerStrs = Vec<KmerStr>;
type WindowKmers = Vec<Kmer>;
type WindowIds = Id;
type WindowTruths = Vec<Truth>;

type BatchKmerStrs = Vec<WindowKmerStrs>;
type BatchKmers = Vec<WindowKmers>;
type BatchIds = Vec<WindowIds>;
type BatchTruths = Vec<WindowTruths>;

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
        let queue = Arc::new(ArrayQueue::new(queue_size));

        // Create KmerWindowGenerator
        let kmer_window_generator =
            KmerWindowGenerator::new(&filename, k, window_size, 0, false, rand);

        let discriminator_masked_generator =
            DiscriminatorMaskedGenerator::new(replacement_pct, k, kmer_window_generator);

        let shutdown = Arc::new(AtomicBool::new(false));
        let exhausted = Arc::new(AtomicBool::new(false));

        let iter = Box::new(discriminator_masked_generator);

        let shutdown_c = Arc::clone(&shutdown);
        let exhausted_c = Arc::clone(&exhausted);

        let queue_handle = Arc::clone(&queue);
        let handle = thread::spawn(move || {
            let mut iter: Box<dyn Iterator<Item = DiscriminatorMasked> + Send> = iter;
            let mut offset = 0;
            let mut rc = false;

            loop {
                let mut finished = false;

                while !finished {

                    let mut window_kmers: Vec<Vec<u8>>;
                    let mut window_id: String;
                    let mut window_truth: Vec<u8>;
    

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

                                let discriminator_masked_generator =
                                    DiscriminatorMaskedGenerator::new(
                                        replacement_pct,
                                        k,
                                        kmer_window_generator,
                                    );

                                iter = Box::new(discriminator_masked_generator);
                                continue;
                            }
                        }
                    };

                    let DiscriminatorMasked { kmers, id, truth } = item;
                    let kmers = kmers
                        .iter()
                        .map(|x| convert_string_to_array(k, x))
                        .collect();
                    window_kmers = kmers;
                    window_id = id;
                    window_truth = truth;

                    let mut batch = (window_kmers, window_id, window_truth);
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

        MaskedKmersGenerator {
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
impl PyIterProtocol for MaskedKmersGenerator {
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
        pyout.set_item("id", result.1).expect("Error with Python");
        pyout
            .set_item("truth", result.2)
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
        let mut batch_kmers: Vec<Vec<Vec<u8>>> = Vec::with_capacity(mypyself.batch_size);
        let mut batch_id: Vec<String> = Vec::with_capacity(mypyself.batch_size);
        let mut batch_truth: Vec<Vec<u8>> = Vec::with_capacity(mypyself.batch_size);

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
                let kmers: Vec<Vec<u8>> = kmers
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
                let kmers: Vec<Vec<u8>> = kmers
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
    m.add_wrapped(wrap_pyfunction!(convert_fasta_to_sfasta))?;
    m.add_wrapped(wrap_pyfunction!(get_headers_from_sfasta))?;
    m.add_wrapped(wrap_pyfunction!(test_sfasta))?;
    m.add_wrapped(wrap_pyfunction!(index_sfasta))?;

    Ok(())
}
