// use pyo3::prelude::*;
// use pyo3::wrap_pyfunction;

// use std::fs::File;
// use std::io::BufReader;

// use finalfusion::prelude::*;


#[cfg(test)]
mod test {
    use std::fs::File;
    use std::io::BufReader;

    use finalfusion::prelude::*;

    #[test]
    fn load_embeddings() {
        let mut reader = BufReader::new(File::open("/home/josephguhlin/work/oracular/test.embed").unwrap());

        // Read the embeddings.
        let embeddings: Embeddings<VocabWrap, StorageWrap> =
            Embeddings::read_embeddings(&mut reader)
            .unwrap();

        assert_eq!(embeddings.dims(), 32, "Test embeddings should have 32 dimensions");

        // println!("{:#?}", embeddings.embedding("ACTGCCCATACGG").unwrap().as_view());

    }
}
