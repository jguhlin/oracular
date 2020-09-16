extern crate rayon;

// NOTE: New naming convention
// Rust-y stuff is "iter" Python is "Generator"

use liboracular::fasta::{parse_ctfasta_target_contexts, parse_fasta_kmers_shuffle};
use liboracular::kmers::{
    DiscriminatorMasked, DiscriminatorMaskedGenerator, Gff3Kmers, Gff3KmersIter, KmerCoordsWindow,
    KmerCoordsWindowIter, KmerWindowGenerator,
};
use liboracular::sfasta;
use liboracular::threads::{
    Sequence, SequenceBatch, SequenceBatchKmers, SequenceKmers, SequenceTargetContexts,
    ThreadCommand,
};

use pyo3::prelude::*;
// use pyo3::wrap_pyfunction;
use pyo3::types::{PyDict, PyList};
use pyo3::PyIterProtocol;

// use pyo3::wrap_pyfunction;

use crossbeam::queue::{ArrayQueue, PopError};
use crossbeam::utils::Backoff;
use pyo3::wrap_pyfunction;
use std::sync::{Arc, RwLock};
use std::thread::JoinHandle;

// use std::fs::File;
// use std::io::BufReader;

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
// Non-batch discriminator masked generator
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
                            mypyself.filename.clone(),
                            mypyself.k,
                            mypyself.window_size,
                            mypyself.offset,
                            mypyself.rc,
                        );

                        let iter = Gff3KmersIter::new(
                            mypyself.gff3filename.clone(),
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
    fn new(k: usize, filename: String, window_size: usize, gff3filename: String) -> Self {
        // Create KmerWindowGenerator
        let kmercoords_window_iter =
            KmerCoordsWindowIter::new(filename.clone(), k, window_size, 0, false);

        let iter = Gff3KmersIter::new(gff3filename.clone(), kmercoords_window_iter, k);
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
        }
    }

    fn types(&self) -> Vec<String> {
        self.types.clone()
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
                            mypyself.filename.clone(),
                            mypyself.k,
                            mypyself.window_size,
                            mypyself.offset,
                            mypyself.rc,
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
    ) -> Self {
        // Create KmerWindowGenerator
        let kmer_window_generator =
            KmerWindowGenerator::new(filename.clone(), k, window_size, 0, false);

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
                            mypyself.filename.clone(),
                            mypyself.k,
                            mypyself.window_size,
                            mypyself.offset,
                            mypyself.rc,
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
            KmerWindowGenerator::new(filename.clone(), k, window_size, 0, false);

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

// Offset & RC Wrapper

// Acc2Tax Discriminator Generator
// TODO: Move to Acc2Tax OR Have Acc2Tax put a wrapper around this fn! (Even
// smarter!)

/*
#[pyclass]
struct DiscriminatorMaskedGeneratorWrapperA2T {
    iter: Box<dyn Iterator<Item = DiscriminatorMasked>>, // + Send>,
    batch_size: usize,
    k: usize,
    offset: usize,
}

#[pyproto]
impl PyIterProtocol for DiscriminatorMaskedGeneratorWrapperA2T {
    fn __iter__(mypyself: PyRefMut<Self>) -> PyResult<PyObject> {
        let gil = Python::acquire_gil();
        let py = gil.python();
        Ok(mypyself.into_py(py))
    }

    fn __next__(mut mypyself: PyRefMut<Self>) -> PyResult<Option<PyObject>> {

        // Generate Batch
        let mut batch_kmers: Vec<Vec<Vec<u8>>> = Vec::with_capacity(mypyself.batch_size);
        let mut batch_id: Vec<String> = Vec::with_capacity(mypyself.batch_size);
        let mut batch_taxons: Vec<Vec<usize>> = Vec::with_capacity(mypyself.batch_size);
        let mut batch_taxon: Vec<usize> = Vec::with_capacity(mypyself.batch_size);
        let mut batch_truth: Vec<Vec<u8>> = Vec::with_capacity(mypyself.batch_size);

        for _ in 0..mypyself.batch_size {
            let item = match mypyself.iter.next() {
                Some(x) => x,
                None    => return Ok(None)
            };

            let DiscriminatorMasked { kmers, id, taxons, taxon, truth} = item;
            let kmers = kmers.iter().map(|x| convert_string_to_array(mypyself.k, x)).collect();
            batch_kmers.push(kmers);
            batch_id.push(id);
            batch_taxons.push(taxons);
            batch_taxon.push(taxon);
            batch_truth.push(truth);
        }

        let gil = Python::acquire_gil();
        let py = gil.python();
        let pyout = PyDict::new(py);
        pyout.set_item("kmers",  batch_kmers ).expect("Error with Python");
        pyout.set_item("id",     batch_id    ).expect("Error with Python");
        pyout.set_item("taxons", batch_taxons).expect("Error with Python");
        pyout.set_item("taxon",  batch_taxon ).expect("Error with Python");
        pyout.set_item("truth",  batch_truth ).expect("Error with Python");
        Ok(Some(pyout.to_object(py)))
    }
}

#[pymethods]
impl DiscriminatorMaskedGeneratorWrapperA2T {
    #[new]
    fn new(
        k: usize,
        filename: String,
        window_size: usize,
        batch_size: usize,
        replacement_pct: f32,
    ) -> Self
    {

        // Create KmerWindowGenerator
        let kmer_window_generator = KmerWindowGenerator::new(
                                        filename,
                                        k.clone(),
                                        window_size,
                                        0);

        let discriminator_masked_generator = DiscriminatorMaskedGenerator::new(
                                        replacement_pct,
                                        k.clone(),
                                        kmer_window_generator);

        DiscriminatorMaskedGeneratorWrapper {
            iter: Box::new(discriminator_masked_generator),
            batch_size: batch_size,
            k,
            offset: 0,
        }
    }
}*/

// SPSC Implementation
// Single producer, single consumer
// One thread processes, the other thread generates the data
// TODO

#[pyclass]
struct CTFasta {
    batch_queue: Arc<ArrayQueue<ThreadCommand<SequenceBatch>>>,
    seq_queue: Arc<ArrayQueue<ThreadCommand<Sequence>>>,
    unshuffled_queue: Arc<ArrayQueue<ThreadCommand<SequenceTargetContexts>>>,
    generator_done: Arc<RwLock<bool>>,
    generator: Option<JoinHandle<()>>,
    batch_size: usize,
    children: Option<Vec<JoinHandle<()>>>,
    seq_queue_shutdown: bool,
    batch_queue_shutdown: bool,
    unshuffled_queue_shutdown: bool,
    num_threads: usize,
}

#[pymethods]
impl CTFasta {
    #[new]
    fn new(
        k: usize,
        filename: String,
        window_size: usize,
        batch_size: usize,
        shuffle_buffer: usize,
        buffer_size: usize,
        num_threads: usize,
    ) -> Self {
        let (
            batch_queue,
            seq_queue,
            unshuffled_queue,
            generator,
            generator_done,
            _children_asleep,
            _jobs,
            children,
        ) = parse_ctfasta_target_contexts(
            k,
            &filename,
            window_size,
            batch_size,
            shuffle_buffer,
            buffer_size,
            num_threads,
        );
        CTFasta {
            batch_queue,
            batch_size,
            seq_queue,
            unshuffled_queue,
            generator: Some(generator),
            generator_done,
            children: Some(children),
            num_threads,
            seq_queue_shutdown: false,
            batch_queue_shutdown: false,
            unshuffled_queue_shutdown: false,
        }
    }

    fn next_batch(&mut self) -> PyResult<PyObject> {
        let mut pop = self.batch_queue.pop();
        let backoff = Backoff::new();
        let batch;

        // Unpark all threads
        for child in self.children.as_ref().unwrap() {
            child.thread().unpark();
        }

        self.generator.as_ref().unwrap().thread().unpark();

        // Generator is done
        if *self.generator_done.read().unwrap() && self.generator.is_some() {
            println!("DBG: Entire file read, shutting down children threads... More data may still be incoming...");
            self.generator
                .take()
                .expect("Unable to get generator thread")
                .join()
                .expect("Unable to join generator thread...");
        }

        // Generator is done...
        if !self.seq_queue_shutdown
            && *self.generator_done.read().unwrap()
            && self.seq_queue.is_empty()
        {
            for _ in 0..self.num_threads {
                match self.seq_queue.push(ThreadCommand::Terminate) {
                    Ok(_) => (),
                    Err(x) => panic!("Unable to send command... {:#?}", x),
                }
                self.seq_queue_shutdown = true;
            }
        }

        if !self.unshuffled_queue_shutdown
            && *self.generator_done.read().unwrap()
            && self.seq_queue.is_empty()
            && self.unshuffled_queue.is_empty()
        //            && self.unshuffled_queue.is_empty()
        {
            match self.unshuffled_queue.push(ThreadCommand::Terminate) {
                Ok(_) => (),
                Err(x) => panic!("Unable to send command... {:#?}", x),
            };
            self.unshuffled_queue_shutdown = true;
        }

        if !self.batch_queue_shutdown
            && *self.generator_done.read().unwrap()
            && self.seq_queue.is_empty()
            && self.unshuffled_queue.is_empty()
            && self.batch_queue.is_empty()
        {
            /* for _ in 0..self.num_threads {
                match self.batch_queue.push(ThreadCommand::Terminate) {
                    Ok(_) => (),
                    Err(x) => panic!("Unable to send command... {:#?}", x)
                }
            } */

            self.batch_queue_shutdown = true;

            for child in self.children.take().expect("Unable to join children") {
                match child.join() {
                    Ok(_) => (),
                    Err(x) => panic!("Error joining worker thread... {:#?}", x),
                }
            }
        }

        if *self.generator_done.read().unwrap()
            && self.seq_queue.is_empty()
            && self.unshuffled_queue.is_empty()
            && self.batch_queue.is_empty()
            && self.batch_queue_shutdown
        {
            println!("Empty!");
            let gil = Python::acquire_gil();
            let py = gil.python();
            let pybatch = PyList::empty(py);
            Ok(pybatch.to_object(py))
        } else {
            while pop == Err(PopError) {
                // Unpark all threads
                for child in self.children.as_ref().unwrap() {
                    child.thread().unpark();
                }

                self.generator.as_ref().unwrap().thread().unpark();

                backoff.snooze();
                pop = self.batch_queue.pop();
            }

            batch = pop.unwrap();

            // let batch = pop.unwrap();

            let gil = Python::acquire_gil();
            let py = gil.python();

            let mut idvec = Vec::with_capacity(self.batch_size);
            let mut targetvec = Vec::with_capacity(self.batch_size);
            let mut contextvec = Vec::with_capacity(self.batch_size);
            // let pybatch = PyList::empty(py);

            for entry in batch.unwrap() {
                // let pyentry = PyDict::new(py);
                idvec.push(entry.id);
                targetvec.push(entry.target);
                contextvec.push(entry.contexts);
                //                pyentry.set_item("ID", entry.id);
                //                pyentry.set_item("target", entry.target);
                // Convert contexts into PyList
                //                let contexts = PyList::new(py,
                // entry.contexts);                
                // pyentry.set_item("contexts", contexts);

                //                batchvec.push(pyentry);
                // pybatch.append(pyentry);
            }

            //Ok(PyList::new(py, batchvec).to_object(py))
            let pyout = PyDict::new(py);
            pyout.set_item("IDs", idvec).expect("Error with Python");
            pyout
                .set_item("contexts", contextvec)
                .expect("Error with Python");
            pyout
                .set_item("targets", targetvec)
                .expect("Error with Python");
            Ok(pyout.to_object(py))
        }
    }
}

#[pyclass]
struct FastaKmersShuffle {
    batch_queue: Arc<ArrayQueue<ThreadCommand<SequenceBatchKmers>>>,
    seq_queue: Arc<ArrayQueue<ThreadCommand<Sequence>>>,
    unshuffled_queue: Arc<ArrayQueue<ThreadCommand<SequenceKmers>>>,
    generator_done: Arc<RwLock<bool>>,
    generator: Option<JoinHandle<()>>,
    batch_size: usize,
    children: Option<Vec<JoinHandle<()>>>,
    seq_queue_shutdown: bool,
    batch_queue_shutdown: bool,
    unshuffled_queue_shutdown: bool,
    num_threads: usize,
}

#[pymethods]
impl FastaKmersShuffle {
    #[new]
    fn new(
        k: usize,
        filename: String,
        sample_size: usize,
        batch_size: usize,
        shuffle_buffer: usize,
        buffer_size: usize,
        num_threads: usize,
    ) -> Self {
        let (batch_queue, seq_queue, unshuffled_queue, generator, generator_done, _jobs, children) =
            parse_fasta_kmers_shuffle(
                k,
                &filename,
                sample_size,
                batch_size,
                shuffle_buffer,
                buffer_size,
                num_threads,
            );
        FastaKmersShuffle {
            batch_queue,
            batch_size,
            seq_queue,
            unshuffled_queue,
            generator: Some(generator),
            generator_done,
            children: Some(children),
            num_threads,
            seq_queue_shutdown: false,
            batch_queue_shutdown: false,
            unshuffled_queue_shutdown: false,
        }
    }

    fn next_batch(&mut self) -> PyResult<PyObject> {
        let mut pop = self.batch_queue.pop();
        let backoff = Backoff::new();
        let batch;

        // Unpark all threads
        for child in self.children.as_ref().unwrap() {
            child.thread().unpark();
        }

        self.generator.as_ref().unwrap().thread().unpark();

        // Generator is done
        if *self.generator_done.read().unwrap() && self.generator.is_some() {
            println!("DBG: Entire file read, shutting down children threads... More data may still be incoming...");
            self.generator
                .take()
                .expect("Unable to get generator thread")
                .join()
                .expect("Unable to join generator thread...");
        }

        // Generator is done...
        if !self.seq_queue_shutdown
            && *self.generator_done.read().unwrap()
            && self.seq_queue.is_empty()
        {
            for _ in 0..self.num_threads {
                match self.seq_queue.push(ThreadCommand::Terminate) {
                    Ok(_) => (),
                    Err(x) => panic!("Unable to send command... {:#?}", x),
                }
                self.seq_queue_shutdown = true;
            }
        }

        if !self.unshuffled_queue_shutdown
            && *self.generator_done.read().unwrap()
            && self.seq_queue.is_empty()
            && self.unshuffled_queue.is_empty()
        //            && self.unshuffled_queue.is_empty()
        {
            match self.unshuffled_queue.push(ThreadCommand::Terminate) {
                Ok(_) => (),
                Err(x) => panic!("Unable to send command... {:#?}", x),
            };
            self.unshuffled_queue_shutdown = true;
        }

        if !self.batch_queue_shutdown
            && *self.generator_done.read().unwrap()
            && self.seq_queue.is_empty()
            && self.unshuffled_queue.is_empty()
            && self.batch_queue.is_empty()
        {
            /* for _ in 0..self.num_threads {
                match self.batch_queue.push(ThreadCommand::Terminate) {
                    Ok(_) => (),
                    Err(x) => panic!("Unable to send command... {:#?}", x)
                }
            } */

            self.batch_queue_shutdown = true;

            for child in self.children.take().expect("Unable to join children") {
                match child.join() {
                    Ok(_) => (),
                    Err(x) => panic!("Error joining worker thread... {:#?}", x),
                }
            }
        }

        if *self.generator_done.read().unwrap()
            && self.seq_queue.is_empty()
            && self.unshuffled_queue.is_empty()
            && self.batch_queue.is_empty()
            && self.batch_queue_shutdown
        {
            println!("Empty!");
            let gil = Python::acquire_gil();
            let py = gil.python();
            let pybatch = PyList::empty(py);
            Ok(pybatch.to_object(py))
        } else {
            while pop == Err(PopError) {
                // Unpark all threads
                for child in self.children.as_ref().unwrap() {
                    child.thread().unpark();
                }

                self.generator.as_ref().unwrap().thread().unpark();

                backoff.snooze();
                pop = self.batch_queue.pop();
            }

            batch = pop.unwrap();

            // let batch = pop.unwrap();

            let gil = Python::acquire_gil();
            let py = gil.python();

            let mut idvec = Vec::with_capacity(self.batch_size);
            let mut kmersvec = Vec::with_capacity(self.batch_size);
            // let pybatch = PyList::empty(py);

            for entry in batch.unwrap() {
                // let pyentry = PyDict::new(py);
                idvec.push(entry.id);
                kmersvec.push(entry.kmers);
                //                pyentry.set_item("ID", entry.id);
                //                pyentry.set_item("target", entry.target);
                // Convert contexts into PyList
                //                let contexts = PyList::new(py,
                // entry.contexts);                
                // pyentry.set_item("contexts", contexts);

                //                batchvec.push(pyentry);
                // pybatch.append(pyentry);
            }

            //Ok(PyList::new(py, batchvec).to_object(py))
            let pyout = PyDict::new(py);
            pyout.set_item("IDs", idvec).expect("Error with Python");
            pyout
                .set_item("kmers", kmersvec)
                .expect("Error with Python");
            Ok(pyout.to_object(py))
        }
    }
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
                            mypyself.filename.clone(),
                            mypyself.k,
                            mypyself.window_size,
                            mypyself.offset,
                            mypyself.rc,
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
        let iter = KmerCoordsWindowIter::new(filename.clone(), k, window_size, 0, start_rc);

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

#[pyfunction]
fn convert_fasta_to_sfasta(input: String, output: String) {
    sfasta::convert_fasta_file(input, output);
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

/// Provides functions for python dealing with Kmers from fasta and sfasta
/// files...
#[pymodule]
fn pyracular(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<CTFasta>()?;
    m.add_class::<FastaKmersGenerator>()?;
    m.add_class::<DiscriminatorMaskedGeneratorWrapper>()?;
    m.add_class::<DiscriminatorMaskedGeneratorWrapperNB>()?;
    m.add_class::<Gff3KmerGenerator>()?;
    m.add_wrapped(wrap_pyfunction!(convert_fasta_to_sfasta))?;
    m.add_wrapped(wrap_pyfunction!(get_headers_from_sfasta))?;
    m.add_wrapped(wrap_pyfunction!(test_sfasta))?;
    m.add_wrapped(wrap_pyfunction!(index_sfasta))?;

    m.add_class::<FastaKmersShuffle>()?; // This one really isn't used anymore, I think

    Ok(())
}
