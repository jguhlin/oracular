use liboracular::fasta::{parse_ntfasta_target_contexts};
use liboracular::threads::{Sequence, ThreadCommand, SequenceBatch, SequenceTargetContexts};

use pyo3::prelude::*;
// use pyo3::wrap_pyfunction;
// use pyo3::{PyIterProtocol, PyClassShell};
use pyo3::types::{PyDict, PyList};

// use pyo3::wrap_pyfunction;

use crossbeam::queue::{ArrayQueue, PopError};
use std::sync::{Arc, RwLock};
use crossbeam::atomic::AtomicCell;
use std::thread::JoinHandle;
use crossbeam::utils::Backoff;


// use std::fs::File;
// use std::io::BufReader;

#[pyclass]
struct NtFasta {
    batch_queue: Arc<ArrayQueue<ThreadCommand<SequenceBatch>>>,
    seq_queue: Arc<ArrayQueue<ThreadCommand<Sequence>>>,
    unshuffled_queue: Arc<ArrayQueue<ThreadCommand<SequenceTargetContexts>>>,
    generator_done: Arc<RwLock<bool>>,
    generator: Option<JoinHandle<()>>,
    jobs: Arc<AtomicCell<usize>>,
    children: Option<Vec<JoinHandle<()>>>,
    seq_queue_shutdown: bool,
    batch_queue_shutdown: bool,
    unshuffled_queue_shutdown: bool,
    num_threads: usize
}

#[pymethods]
impl NtFasta {
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
        let (batch_queue, seq_queue, unshuffled_queue, generator, generator_done, jobs, children) = parse_ntfasta_target_contexts(k, &filename, window_size, batch_size, shuffle_buffer, buffer_size, num_threads);
        NtFasta { 
            batch_queue, 
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

                backoff.snooze();
                backoff.snooze();
                backoff.snooze();
                pop = self.batch_queue.pop();
            }

            batch = pop.unwrap();

            // let batch = pop.unwrap();

            let gil = Python::acquire_gil();
            let py = gil.python();
            let pybatch = PyList::empty(py);

            for entry in batch.unwrap() {
                let pyentry = PyDict::new(py);
                pyentry.set_item("ID", entry.id);
                pyentry.set_item("target", entry.target);
                // Convert contexts into PyList
                let contexts = PyList::new(py, entry.contexts);
                pyentry.set_item("contexts", contexts);

                pybatch.append(pyentry);
            }

            Ok(pybatch.to_object(py))
        }
    }

}


/// This module is a python module implemented in Rust.
#[pymodule]
fn pyracular(py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<NtFasta>()?;
    Ok(())
}


/*
#[cfg(test)]
mod test {
    // use std::fs::File;
    // use std::io::BufReader;

}
*/