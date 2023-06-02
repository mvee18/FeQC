use clap::Parser;
use itertools::Itertools;
use rayon::prelude::*;
use std::collections::HashMap;
use std::error::Error;
use std::fs::read_dir;
use std::fs::File;
use std::fs::OpenOptions;
use std::io::BufRead;
use std::io::Read;
use std::io::{BufReader, Write};
use std::path::PathBuf;

#[derive(Clone, Debug)]
pub struct FastqFile {
    pub fp: String,
    pub records: Vec<Record>,
}

#[derive(Clone, Debug)]
pub struct Record {
    pub id: String,
    pub seq: String,
    pub plus: String,
    pub qual: String,
}

impl Record {
    pub fn new() -> Record {
        Record {
            id: String::new(),
            seq: String::new(),
            plus: String::new(),
            qual: String::new(),
        }
    }
    pub fn verify_integrity(&self) -> bool {
        // Check that the length of the sequence and quality are the same.
        self.seq.len() == self.qual.len()
    }
    pub fn write_record(&self, fp: &PathBuf) -> Result<(), Box<dyn Error>> {
        // Write the record to a file.
        let mut file = OpenOptions::new()
            .create(true)
            .write(true)
            .append(true)
            .open(fp)
            .unwrap();
        writeln!(file, "{}", self.id)?;
        writeln!(file, "{}", self.seq)?;
        writeln!(file, "{}", self.plus)?;
        writeln!(file, "{}", self.qual)?;

        Ok(())
    }
}

pub fn read_fastq_file(path: &str) -> Result<Vec<Record>, Box<dyn Error>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut records = Vec::new();
    let lines = reader.lines().into_iter();
    for (i, item) in lines.enumerate() {
        let line = item?;
        match i % 4 {
            0 => {
                let mut record = Record::new();
                record.id = line;
                records.push(record);
            }
            1 => {
                records.last_mut().unwrap().seq = line;
            }
            2 => {
                records.last_mut().unwrap().plus = line;
            }
            3 => {
                records.last_mut().unwrap().qual = line;
                // Verify that the record is valid.
                if !records.last().unwrap().verify_integrity() {
                    println!("Invalid record on line {}", i);
                    return Err("Invalid record".into());
                }
            }
            _ => {}
        }
    }
    Ok(records)
}

pub fn convert_ascii_to_score(ascii: char) -> f64 {
    ascii as u8 as f64 - 33.0
}

pub fn calculate_fastq_quality_score(qual: &str) -> f64 {
    let mut sum = 0.0;
    for c in qual.chars() {
        sum += convert_ascii_to_score(c)
    }
    sum / qual.len() as f64
}

pub fn get_average_quality_score(records: &Vec<Record>) -> f64 {
    let mut sum = 0.0;
    for record in records {
        sum += calculate_fastq_quality_score(&record.qual);
    }
    sum / records.len() as f64
}

pub fn rayon_get_average_quality_score(records: &Vec<Record>) -> f64 {
    let average = records
        .par_iter()
        .map(|r| calculate_fastq_quality_score(&r.qual))
        .sum::<f64>()
        / records.len() as f64;

    average
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input directory path.
    #[arg(short, long)]
    input: String,
}

pub fn read_and_get_records() -> Result<Vec<FastqFile>, Box<dyn Error>> {
    // Add clap argument to read file path.
    let args: Args = Args::parse();

    // We want to read all the files in the directory.
    let paths = read_dir(&args.input).unwrap();

    // Initialize fastqfile vec.
    let mut fastqfiles: Vec<FastqFile> = Vec::new();

    for p in paths {
        let path = p.unwrap().path();
        let path_str = path.to_str().unwrap();
        if path_str.ends_with(".fastq") {
            let records = match read_fastq_file(path_str) {
                Ok(r) => r,
                Err(e) => {
                    println!("Error reading file {}: {}", path_str, e);
                    return Err(e);
                }
            };

            let fastqfile = FastqFile {
                fp: path_str.to_string(),
                records: records,
            };

            fastqfiles.push(fastqfile);
        }
    }

    if fastqfiles.len() == 0 {
        return Err("No fastq files found".into());
    } else {
        Ok(fastqfiles)
    }
}

fn get_index_scores_and_count(records: &Vec<Record>) -> (HashMap<i32, i32>, HashMap<i32, f64>) {
    let mut index_count: HashMap<i32, i32> = HashMap::new();

    // What if the first sequence isn't the longest?
    let mut index_qual: HashMap<i32, f64> = HashMap::new();

    for record in records {
        println!("{}", record.qual.len());
        for (i, c) in record.qual.chars().enumerate() {
            index_qual
                .entry(i as i32)
                .and_modify(|e| *e += convert_ascii_to_score(c))
                .or_insert(convert_ascii_to_score(c));

            index_count
                .entry(i as i32)
                .and_modify(|e| *e += 1)
                .or_insert(1);
        }
    }

    return (index_count, index_qual);
}

/// This will calculate the average quality score across all indices.
/// # Arguments
/// * `records` - A vector of records.
/// # Returns
/// * A vector of f64s representing the average quality score at each index.
pub fn average_quality_at_index(records: &Vec<Record>) -> Vec<f64> {
    let (index_count, index_qual) = get_index_scores_and_count(records);

    let mut result: Vec<f64> = vec![0.0; index_count.len()];

    for k in index_qual.keys() {
        result[*k as usize] = index_qual[k] / index_count[k] as f64;
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mock_read_and_get_records() -> Result<Vec<FastqFile>, Box<dyn Error>> {
        let f1_fp = "src/test_fastq/Staph_R1_head.fastq";
        let f2_fp = "src/test_fastq/Staph_R2_head.fastq";

        let fq1 = FastqFile {
            fp: f1_fp.to_string(),
            records: read_fastq_file(f1_fp).unwrap(),
        };

        let fq2 = FastqFile {
            fp: f2_fp.to_string(),
            records: read_fastq_file(f2_fp).unwrap(),
        };

        Ok(vec![fq1, fq2])
    }

    fn hashmap_comparer(
        a: &HashMap<i32, i32>,
        b: &HashMap<i32, i32>,
    ) -> Result<bool, Box<dyn Error>> {
        if a.len() != b.len() {
            return Err(format!("{} != {} not equal", a.len(), b.len()).into());
        }

        for k in a.keys() {
            if a[k] != b[k] {
                return Err("Values not equal".into());
            }
        }

        Ok(true)
    }

    #[test]
    fn test_index_count() {
        let fqs = mock_read_and_get_records().unwrap();

        // The wanted index count has keys 1-62 with value 2 and key values 62-151 with value 1.
        let wanted_index_count_f1: HashMap<i32, i32> = (0..63)
            .map(|i| (i, 2))
            .chain((63..150).map(|i| (i, 1)))
            .collect();

        let wanted_index_count_f2: HashMap<i32, i32> = (0..68)
            .map(|i| (i, 2))
            .chain((68..151).map(|i| (i, 1)))
            .collect();

        for (c, fq) in fqs.iter().enumerate() {
            println!("{}", c);
            println!("{:?}", fq.fp);
            let (index_count, _) = get_index_scores_and_count(&fq.records);
            if c == 0 {
                assert!(hashmap_comparer(&index_count, &wanted_index_count_f1).unwrap());
            } else {
                assert!(hashmap_comparer(&index_count, &wanted_index_count_f2).unwrap());
            }
        }
    }

    #[test]
    fn test_average_index_qual() {
        // We can read from the r1 and r2 wanted files for what our vector should be.
        let mut f1 = File::open("src/test_fastq/r1_wanted_avg_index.txt").unwrap();
        let mut f2 = File::open("src/test_fastq/r2_wanted_avg_index.txt").unwrap();

        // The file is line by line with a single int value.
        let mut wanted_f1: Vec<f64> = Vec::new();
        let mut wanted_f2: Vec<f64> = Vec::new();

        let mut buf = String::new();
        f1.read_to_string(&mut buf).unwrap();
        for line in buf.lines() {
            wanted_f1.push(line.parse::<f64>().unwrap());
        }

        buf.clear();

        f2.read_to_string(&mut buf).unwrap();
        for line in buf.lines() {
            wanted_f2.push(line.parse::<f64>().unwrap());
        }

        let fqs = mock_read_and_get_records().unwrap();
        let f1_avg_index_qual = average_quality_at_index(&fqs[0].records);
        let f2_avg_index_qual = average_quality_at_index(&fqs[1].records);

        assert_eq!(f1_avg_index_qual.len(), wanted_f1.len());
        assert_eq!(f2_avg_index_qual.len(), wanted_f2.len());
    }
}
