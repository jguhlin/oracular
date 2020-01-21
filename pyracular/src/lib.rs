use liboracular::fasta::{parse_ntfasta_target_contexts};
use liboracular::threads::{ThreadCommand, SequenceBatch};

use pyo3::prelude::*;
// use pyo3::wrap_pyfunction;
// use pyo3::{PyIterProtocol, PyClassShell};
use pyo3::types::{PyDict, PyList};

use pyo3::wrap_pyfunction;

use crossbeam::queue::{ArrayQueue, PushError};
use std::sync::{Arc, RwLock};
use crossbeam::atomic::AtomicCell;
use std::thread::JoinHandle;


// use std::fs::File;
// use std::io::BufReader;

#[pyclass]
struct NtFasta {
    batch_queue: Arc<ArrayQueue<ThreadCommand<SequenceBatch>>>,
    generator_done: Arc<RwLock<bool>>,
    jobs: Arc<AtomicCell<usize>>,
    children: Vec<JoinHandle<()>>,
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
        let (batch_queue, generator_done, jobs, children) = parse_ntfasta_target_contexts(k, &filename, window_size, batch_size, shuffle_buffer, buffer_size, num_threads);
        NtFasta { batch_queue, generator_done, jobs, children }
    }

    fn next_batch(&self) -> PyResult<PyObject> {
        let batch = self.batch_queue.pop().unwrap(); // TODO: Add code for when we run out of sequence

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