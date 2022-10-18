use clap::Parser;
use std::error::Error;
use std::fs::read_dir;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::path::Path;

pub struct Record {
    id: String,
    seq: String,
    plus: String,
    qual: String,
}

impl Record {
    fn new() -> Record {
        Record {
            id: String::new(),
            seq: String::new(),
            plus: String::new(),
            qual: String::new(),
        }
    }
    fn verify_integrity(&self) -> bool {
        // Check that the length of the sequence and quality are the same.
        self.seq.len() == self.qual.len()
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

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input directory path.
    #[arg(short, long)]
    input: String,
}

fn main() {
    // Add clap argument to read file path.
    let args: Args = Args::parse();

    // We want to read all the files in the directory.
    let paths = read_dir(&args.input).unwrap();
    for p in paths {
        let path = p.unwrap().path();
        let path_str = path.to_str().unwrap();
        if path_str.ends_with(".fastq") {
            let records = match read_fastq_file(path_str) {
                Ok(r) => r,
                Err(e) => {
                    println!("Error reading file {}: {}", path_str, e);
                    return;
                }
            };
            let avg = get_average_quality_score(&records);
            let fp = Path::new(path_str);

            println!("{}\t{}", fp.file_name().unwrap().to_str().unwrap(), avg);
        }
    }
}
