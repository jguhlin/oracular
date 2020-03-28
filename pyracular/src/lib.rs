extern crate rayon;

use liboracular::kmers::{KmerWindowGenerator, DiscriminatorMasked, DiscriminatorMaskedGenerator};
use liboracular::fasta::{parse_ctfasta_target_contexts, parse_fasta_kmers};
use liboracular::threads::{Sequence, ThreadCommand, SequenceBatch, SequenceTargetContexts, SequenceBatchKmers, SequenceKmers};
use liboracular::sfasta;

use pyo3::prelude::*;
// use pyo3::wrap_pyfunction;
use pyo3::types::{PyDict, PyList};
use pyo3::{PyIterProtocol};

// use pyo3::wrap_pyfunction;

use crossbeam::queue::{ArrayQueue, PopError};
use std::sync::{Arc, RwLock};
use crossbeam::atomic::AtomicCell;
use std::thread::JoinHandle;
use crossbeam::utils::Backoff;
use pyo3::wrap_pyfunction;


// use std::fs::File;
// use std::io::BufReader;

#[inline(always)]
fn convert_string_to_array(k: usize, s: &[u8]) -> Vec<u8> {
    let mut out: Vec<u8> = vec![0; k*5];

    let mut x = 0;
    for c in s {
        match c {
            65 => out[5*x] = 1,   // A
            84 => out[5*x+1] = 1, // T
            67 => out[5*x+3] = 1, // C
            71 => out[5*x+4] = 1, // G
            78 => out[5*x+2] = 1, // N
            _  => out[5*x+2] = 1, // N for everything else...
            // A_  => { out[5*x+2] = 1; println!("Invalid Character! {} in {}", c, std::str::from_utf8(s).unwrap()) }  // N
        };

        x += 1;
    }

    out
}

#[pyclass]
struct DiscriminatorMaskedGeneratorWrapper {
    iter: Box<dyn Iterator<Item = DiscriminatorMasked>>, // + Send>,
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
                None    => { 
                    mypyself.offset = mypyself.offset + 1;
                    if (mypyself.k == mypyself.offset) && mypyself.rc {
                        return Ok(None)
                    } else {

                        if mypyself.k == mypyself.offset {
                            mypyself.rc = true;
                            mypyself.k = 0;
                        }

                        let kmer_window_generator = KmerWindowGenerator::new(
                            mypyself.filename.clone(), 
                            mypyself.k.clone(), 
                            mypyself.window_size.clone(),
                            mypyself.offset.clone(),
                            mypyself.rc.clone(),
                        );

                        let discriminator_masked_generator = DiscriminatorMaskedGenerator::new(
                            mypyself.replacement_pct.clone(),
                            mypyself.k.clone(),
                            kmer_window_generator);

                        mypyself.iter = Box::new(discriminator_masked_generator);
                        continue
                    }
                }
            };

            let DiscriminatorMasked { kmers, id, truth} = item;
            let kmers = kmers.iter().map(|x| convert_string_to_array(mypyself.k, x)).collect();
            batch_kmers.push(kmers);
            batch_id.push(id);
            batch_truth.push(truth);
        }

        let gil = Python::acquire_gil();
        let py = gil.python();
        let pyout = PyDict::new(py);
        pyout.set_item("kmers",  batch_kmers ).expect("Error with Python");
        pyout.set_item("id",     batch_id    ).expect("Error with Python");
        pyout.set_item("truth",  batch_truth ).expect("Error with Python");
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
    ) -> Self 
    {

        // Create KmerWindowGenerator
        let kmer_window_generator = KmerWindowGenerator::new(
                                        filename.clone(), 
                                        k.clone(), 
                                        window_size,
                                        0,
                                        false);
        
        let discriminator_masked_generator = DiscriminatorMaskedGenerator::new(
                                        replacement_pct,
                                        k.clone(),
                                        kmer_window_generator);

        DiscriminatorMaskedGeneratorWrapper { 
            iter: Box::new(discriminator_masked_generator),
            batch_size: batch_size,
            k,
            offset: 0,
            filename: filename.clone(),
            window_size: window_size.clone(),
            replacement_pct,
            rc: false,
        }
    }
}

// Acc2Tax Discriminator Generator
// TODO: Move to Acc2Tax OR Have Acc2Tax put a wrapper around this fn! (Even smarter!)

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
    jobs: Arc<AtomicCell<usize>>,
    batch_size: usize,
    children_asleep: Arc<RwLock<bool>>,
    children: Option<Vec<JoinHandle<()>>>,
    seq_queue_shutdown: bool,
    batch_queue_shutdown: bool,
    unshuffled_queue_shutdown: bool,
    num_threads: usize
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
        num_threads: usize) -> Self 
    {
        let (batch_queue, seq_queue, unshuffled_queue, generator, generator_done, children_asleep, jobs, children) = parse_ctfasta_target_contexts(k, &filename, window_size, batch_size, shuffle_buffer, buffer_size, num_threads);
        CTFasta { 
            batch_queue, 
            batch_size,
            seq_queue, 
            unshuffled_queue,
            generator: Some(generator), 
            generator_done, 
            children_asleep,
            jobs, 
            children: Some(children), 
            num_threads,
            seq_queue_shutdown: false,
            batch_queue_shutdown: false,
            unshuffled_queue_shutdown: false
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
        if *self.generator_done.read().unwrap() {
            if self.generator.is_some() {
                println!("DBG: Entire file read, shutting down children threads... More data may still be incoming...");
                self.generator.take().expect("Unable to get generator thread").join().expect("Unable to join generator thread...");
            }
        }

        // Generator is done...
        if  !self.seq_queue_shutdown 
            && *self.generator_done.read().unwrap() 
            && self.seq_queue.is_empty() 
            
        {
            for _ in 0..self.num_threads {
                match self.seq_queue.push(ThreadCommand::Terminate) {
                    Ok(_) => (),
                    Err(x) => panic!("Unable to send command... {:#?}", x)
                }
            self.seq_queue_shutdown = true;
            }
        }

        if  !self.unshuffled_queue_shutdown 
            && *self.generator_done.read().unwrap() 
            && self.seq_queue.is_empty()
            && self.unshuffled_queue.is_empty()
//            && self.unshuffled_queue.is_empty() 
            
        {
            match self.unshuffled_queue.push(ThreadCommand::Terminate) {
                Ok(_) => (),
                Err(x) => panic!("Unable to send command... {:#?}", x)
            };
            self.unshuffled_queue_shutdown = true;
        }

        if  !self.batch_queue_shutdown 
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
                    Err(x) => panic!("Error joining worker thread... {:#?}", x)
                }
            }
        }

        if  *self.generator_done.read().unwrap() 
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
//                let contexts = PyList::new(py, entry.contexts);
//                pyentry.set_item("contexts", contexts);

//                batchvec.push(pyentry);
                // pybatch.append(pyentry);
            }

            //Ok(PyList::new(py, batchvec).to_object(py))
            let pyout = PyDict::new(py);
            pyout.set_item("IDs", idvec).expect("Error with Python");
            pyout.set_item("contexts", contextvec).expect("Error with Python");
            pyout.set_item("targets", targetvec).expect("Error with Python");
            Ok(pyout.to_object(py))
        }
    }

}

#[pyclass]
struct FastaKmers {
    batch_queue: Arc<ArrayQueue<ThreadCommand<SequenceBatchKmers>>>,
    seq_queue: Arc<ArrayQueue<ThreadCommand<Sequence>>>,
    unshuffled_queue: Arc<ArrayQueue<ThreadCommand<SequenceKmers>>>,
    generator_done: Arc<RwLock<bool>>,
    generator: Option<JoinHandle<()>>,
    jobs: Arc<AtomicCell<usize>>,
    batch_size: usize,
    children: Option<Vec<JoinHandle<()>>>,
    seq_queue_shutdown: bool,
    batch_queue_shutdown: bool,
    unshuffled_queue_shutdown: bool,
    num_threads: usize
}

#[pymethods]
impl FastaKmers {
    #[new]
    fn new(
        k: usize, 
        filename: String, 
        sample_size: usize, 
        batch_size: usize, 
        shuffle_buffer: usize, 
        buffer_size: usize,
        num_threads: usize) -> Self 
    {
        let (batch_queue, seq_queue, unshuffled_queue, generator, generator_done, jobs, children) = parse_fasta_kmers(k, &filename, sample_size, batch_size, shuffle_buffer, buffer_size, num_threads);
        FastaKmers { 
            batch_queue, 
            batch_size,
            seq_queue, 
            unshuffled_queue,
            generator: Some(generator), 
            generator_done, 
            jobs, 
            children: Some(children), 
            num_threads,
            seq_queue_shutdown: false,
            batch_queue_shutdown: false,
            unshuffled_queue_shutdown: false
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
        if *self.generator_done.read().unwrap() {
            if self.generator.is_some() {
                println!("DBG: Entire file read, shutting down children threads... More data may still be incoming...");
                self.generator.take().expect("Unable to get generator thread").join().expect("Unable to join generator thread...");
            }
        }

        // Generator is done...
        if  !self.seq_queue_shutdown 
            && *self.generator_done.read().unwrap() 
            && self.seq_queue.is_empty() 
            
        {
            for _ in 0..self.num_threads {
                match self.seq_queue.push(ThreadCommand::Terminate) {
                    Ok(_) => (),
                    Err(x) => panic!("Unable to send command... {:#?}", x)
                }
            self.seq_queue_shutdown = true;
            }
        }

        if  !self.unshuffled_queue_shutdown 
            && *self.generator_done.read().unwrap() 
            && self.seq_queue.is_empty()
            && self.unshuffled_queue.is_empty()
//            && self.unshuffled_queue.is_empty() 
            
        {
            match self.unshuffled_queue.push(ThreadCommand::Terminate) {
                Ok(_) => (),
                Err(x) => panic!("Unable to send command... {:#?}", x)
            };
            self.unshuffled_queue_shutdown = true;
        }

        if  !self.batch_queue_shutdown 
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
                    Err(x) => panic!("Error joining worker thread... {:#?}", x)
                }
            }
        }

        if  *self.generator_done.read().unwrap() 
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
//                let contexts = PyList::new(py, entry.contexts);
//                pyentry.set_item("contexts", contexts);

//                batchvec.push(pyentry);
                // pybatch.append(pyentry);
            }

            //Ok(PyList::new(py, batchvec).to_object(py))
            let pyout = PyDict::new(py);
            pyout.set_item("IDs", idvec).expect("Error with Python");
            pyout.set_item("kmers", kmersvec).expect("Error with Python");
            Ok(pyout.to_object(py))
        }
    }

}

#[pyfunction]
fn convert_fasta_to_sfasta(input: String, output: String) {
    sfasta::convert_fasta_file(input, output);
}

#[pyfunction]
fn get_headers_from_sfasta(input: String) -> Vec<String>
{
    sfasta::get_headers_from_sfasta(input)
}

/// This module is a python module implemented in Rust.
#[pymodule]
fn pyracular(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<CTFasta>()?;
    m.add_class::<FastaKmers>()?;
    m.add_class::<DiscriminatorMaskedGeneratorWrapper>()?;
    m.add_wrapped(wrap_pyfunction!(convert_fasta_to_sfasta))?;
    m.add_wrapped(wrap_pyfunction!(get_headers_from_sfasta))?;

    Ok(())
}