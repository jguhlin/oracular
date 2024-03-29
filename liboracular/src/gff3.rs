use hashbrown::HashMap;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;

use crate::intervals;

use linked_hash_set::LinkedHashSet;
use rust_lapper::{Interval, Lapper};

use itertools::Itertools;

#[derive(Debug, PartialEq)]
pub enum Strand {
    Plus,  // +
    Minus, // -
}

#[derive(Debug)]
pub struct Gff3Entry {
    pub landmark: String,
    pub source: String,
    pub feature_type: String,
    pub feature_type_strand: String,
    pub start: usize,
    pub end: usize,
    pub score: Option<f32>,
    pub strand: Option<Strand>,
    pub phase: Option<u8>,  //None, 0, 1, or 2
    pub attributes: String, // HashMap<String, String>, // TODO: Implement this...
}

fn is_blank(x: &str) -> bool {
    x == "." || x.is_empty()
}

fn parse_score(x: &str) -> Option<f32> {
    if is_blank(x) {
        return None;
    }
    Some(x.parse::<f32>().expect("Unable to parse score:"))
}

fn parse_strand(x: &str) -> Option<Strand> {
    if is_blank(x) {
        return None;
    }
    if x == "+" {
        Some(Strand::Plus)
    } else if x == "-" {
        Some(Strand::Minus)
    } else {
        None
    }
}

fn parse_phase(x: &str) -> Option<u8> {
    if is_blank(x) {
        return None;
    }
    let n: u8 = x.parse().expect("Unable to parse phase");
    if n > 2 {
        None
    } else {
        Some(n)
    }
}

pub fn parse_gff3_file(filename: &str) -> (Vec<Gff3Entry>, LinkedHashSet<String>) {
    let file = File::open(filename).expect("Unable to open file");
    let reader = BufReader::new(file);
    let lines = reader.lines();

    let mut entries = Vec::with_capacity(1024 * 1024);
    let mut types = LinkedHashSet::new();

    for line in lines {
        if let Some(x) = parse_gff3_line(&line.expect("Unable to parse line")) {
            types.insert(x.feature_type_strand.clone());
            entries.push(x)
        }
    }

    (entries, types)
}

pub fn start_stop_sort(start: usize, stop: usize) -> (usize, usize) {
    let mut v = vec![start, stop];
    v.sort_unstable();
    (v[0], v[1])
}

pub fn get_gff3_intervals(filename: &str) -> (intervals::IntervalMap<Vec<u8>>, Vec<String>) {
    let (entries, types) = parse_gff3_file(filename);

    let types: Vec<String> = types.into_iter().collect();

    let mut types_one_hot: HashMap<String, Vec<u8>> = Default::default();

    let types_len = types.len(); // Calculate only once...

    // Create a one-hot map for each type found in the GFF file
    for (x, t) in types.iter().enumerate() {
        let mut b: Vec<u8> = vec![0; types_len];
        b[x] = 1;

        types_one_hot.insert(t.clone(), b);
    }

    let entries = entries
        .into_iter()
        .map(|x| (x.landmark.clone(), x)) // Convert to tuples of (landmark, x)
        .into_group_map(); // Put it into a hashmap of landmark -> Vec(xs)

    let mut intervals: intervals::IntervalMap<Vec<u8>> = Default::default();

    for (landmark, vals) in entries.iter() {
        let data = vals.iter().map(|x| Interval {
            start: x.start as u32,
            stop: x.end as u32,
            val: types_one_hot
                .get(&x.feature_type_strand)
                .expect("Missing type!")
                .clone(),
        });
        let mut lapper = Lapper::new(data.collect());
        lapper.merge_overlaps();
        intervals.landmarks.insert(landmark.to_string(), lapper);
    }

    (intervals, types)
}

pub fn parse_gff3_line(line: &str) -> Option<Gff3Entry> {
    if line.starts_with('#') {
        // Line is a comment...
        return None;
    }
    let x: Vec<&str> = line.split('\t').collect();
    assert!(
        x.len() == 9,
        "GFF3 has invalid entity count: found {} expected 9",
        x.len()
    );

    let (start, end) = start_stop_sort(
        x[3].parse().expect("Unable to parse start for entry"),
        x[4].parse().expect("Unable to parse end for entry"),
    );

    let feature_type = x[2].to_string().to_lowercase();
    // Filter out regions (Such as 1..end being a "chromosome")
    // to_lowercase because formatting for GFF3 files is iffy
    if feature_type == "region" {
        return None;
    }
    if feature_type == "chromosome" {
        return None;
    }
    if feature_type == "scaffold" {
        return None;
    }
    if feature_type == "contig" {
        return None;
    }
    // Probably will need more...

    let strand = parse_strand(x[6]);
    let feature_type_strand = match strand {
        Some(Strand::Plus) => format!("{}_Plus", feature_type),
        Some(Strand::Minus) => format!("{}_Minus", feature_type),
        None => feature_type.clone(),
    };

    Some(Gff3Entry {
        landmark: x[0].to_string(),
        source: x[1].to_string(),
        feature_type,
        feature_type_strand,
        start,
        end,
        score: parse_score(x[5]),
        strand,
        phase: parse_phase(x[7]),
        attributes: x[8].to_string(),
    })
}

// TODO: Add a test, make sure GFF3 is the same as the intervals we find...
#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;

    #[test]
    pub fn test_parse_score() {
        assert!(parse_score("") == None);
        assert!(parse_score(".") == None);
        assert!(parse_score("1.0") == Some(1.0));
        assert!(parse_score("1") == Some(1.0));
    }

    #[test]
    pub fn test_parse_strand() {
        assert!(parse_strand("+") == Some(Strand::Plus));
        assert!(parse_strand("-") == Some(Strand::Minus));
        assert!(parse_strand("") == None);
        assert!(parse_strand("Notastrand") == None);
    }

    #[test]
    pub fn test_parse_phase() {
        assert!(parse_phase("") == None);
        assert!(parse_phase("0") == Some(0));
        assert!(parse_phase("1") == Some(1));
        assert!(parse_phase("2") == Some(2));
        assert!(parse_phase("3") == None);
    }

    #[test]
    pub fn test_gff3_parse() {
        let mut file = File::open("test_data/Dmel_head30.gff3")
            .expect("Unable to open file test_data/dmel.gff3");
        let mut contents = String::new();
        file.read_to_string(&mut contents)
            .expect("Unable to read file into string");

        for line in contents.lines() {
            let x = parse_gff3_line(&line);
            println!("{:#?}", x);
        }
    }

    #[test]
    #[should_panic]
    pub fn test_gff3_too_many_entries() {
        parse_gff3_line(
            "NC_004354.4	RefSeq	exon	124370	125409	.	-	.	ID=exon .   .   .   .   .   .\t.\t.\t",
        );
    }

    #[test]
    #[should_panic]
    pub fn test_gff3_not_enough_entries() {
        parse_gff3_line("NC_004354.4	RefSeq	exon	124370	125409");
    }

    #[test]
    pub fn test_gff3_nones() {
        assert!(parse_gff3_line(
            "NC_004354.4	RefSeq	ChRoMOsOme	124370	125409	.	-	.	ID=exon .   .   .   .   .   ."
        )
        .is_none());
        assert!(parse_gff3_line(
            "NC_004354.4	RefSeq	chromosome	124370	125409	.	-	.	ID=exon .   .   .   .   .   ."
        )
        .is_none());
        assert!(parse_gff3_line(
            "NC_004354.4	RefSeq	CHROMOSOME	124370	125409	.	-	.	ID=exon .   .   .   .   .   ."
        )
        .is_none());
        assert!(parse_gff3_line(
            "NC_004354.4	RefSeq	Chromosome	124370	125409	.	-	.	ID=exon .   .   .   .   .   ."
        )
        .is_none());
        assert!(parse_gff3_line(
            "NC_004354.4	RefSeq	Region	124370	125409	.	-	.	ID=exon .   .   .   .   .   ."
        )
        .is_none());
        assert!(parse_gff3_line(
            "NC_004354.4	RefSeq	scaffold	124370	125409	.	-	.	ID=exon .   .   .   .   .   ."
        )
        .is_none());
        assert!(parse_gff3_line(
            "NC_004354.4	RefSeq	contig	124370	125409	.	-	.	ID=exon .   .   .   .   .   ."
        )
        .is_none());
    }

    #[test]
    pub fn test_gff3_feature_types() {
        let x =
            parse_gff3_line("NC_004354.4	RefSeq	exon	124370	125409	.	-	.	ID=exon .   .   .   .   .   .");
        assert!(x.unwrap().feature_type_strand.ends_with("_Minus"));

        let x =
            parse_gff3_line("NC_004354.4	RefSeq	exon	124370	125409	.	+	.	ID=exon .   .   .   .   .   .");
        assert!(x.unwrap().feature_type_strand.ends_with("_Plus"));

        let x =
            parse_gff3_line("NC_004354.4	RefSeq	exon	124370	125409	.	.	.	ID=exon .   .   .   .   .   .");
        assert!(x.unwrap().feature_type_strand == "exon");
    }

    #[test]
    pub fn test_strand_impl() {
        assert!(Strand::Minus != Strand::Plus);
    }

    #[test]
    pub fn test_gff3_to_intervals() {
        let (intervals, _) = get_gff3_intervals("test_data/Dmel_head30.gff3");

        let first_landmark = intervals.landmarks.keys().next().unwrap();

        println!("{:#?}", intervals.landmarks.get(first_landmark).unwrap());

        assert!(first_landmark == "NC_004354.4");
        println!(
            "{:#?}",
            intervals.landmarks.get(first_landmark).unwrap().cov()
        );
        assert!(21668 == intervals.landmarks.get(first_landmark).unwrap().cov());
    }
}
