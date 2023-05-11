use clap::Parser;
use rayon::prelude::*;
use std::error::Error;
use std::fs::read_dir;
use std::fs::File;
use std::fs::OpenOptions;
use std::io::BufRead;
use std::io::{BufReader, Write};
use std::path::PathBuf;

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

pub fn calculate_fastq_quality_score(qual: &str) -> f64 {
    let mut sum = 0.0;
    for c in qual.chars() {
        sum += c as u8 as f64 - 33.0;
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
