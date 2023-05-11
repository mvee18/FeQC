// This is a test of Rust's fan in, fan out patterns.
#![allow(unused_imports)]
use fefastq::{calculate_fastq_quality_score, Record};
use rand::Rng;
use rayon::prelude::*;
use std::error::Error;
use std::io;
use std::path::Path;
use std::sync::mpsc::{self, channel, Receiver};
use std::thread;
use std::time::{Duration, Instant};
use threadpool::ThreadPool;

#[allow(dead_code)]
fn basic_thread() {
    let test_vec: Vec<&str> = vec!["a", "b", "c", "d", "e", "f", "g", "h"];
    let handle = thread::spawn(move || {
        for i in 1..3 {
            println!("Hello from the {} thread!", i);
            println!("Here's the vector: {:?}", test_vec);
        }
    });

    handle.join().unwrap();
}

#[allow(dead_code)]
fn channeling() {
    let (tx, rx) = mpsc::channel();

    let tx1 = tx.clone();
    thread::spawn(move || {
        let vals = vec![
            String::from("hi"),
            String::from("from"),
            String::from("the"),
            String::from("thread"),
        ];

        for val in vals {
            tx1.send(val).unwrap();
            thread::sleep(Duration::from_secs(1));
        }
    });

    thread::spawn(move || {
        let vals = vec![
            String::from("more"),
            String::from("messages"),
            String::from("for"),
            String::from("you"),
        ];

        for val in vals {
            tx.send(val).unwrap();
            thread::sleep(Duration::from_secs(1));
        }
    });

    for received in rx {
        println!("Got: {}", received);
    }
}

fn one_to_one_mpsc(records: Vec<Record>) -> f64 {
    let start = Instant::now();
    let (tx, rx) = mpsc::channel();
    let mut sum = 0.0;
    let mut count = 0;

    thread::spawn(move || {
        println!("Hello from the thread!");
        for record in records {
            let val = calculate_fastq_quality_score(&record.qual);
            tx.send(val).unwrap();
        }
    });

    for rec in rx {
        sum += rec;
        count += 1;
    }

    let average = sum / count as f64;

    let duration = start.elapsed();
    println!("Time elapsed in one_to_one() is: {:?}", duration);
    println!("Average quality score: {}", average);
    average
}

fn sequential(records: &Vec<Record>) -> f64 {
    let start = Instant::now();

    let average = records
        .iter()
        .map(|r| calculate_fastq_quality_score(&r.qual))
        .sum::<f64>()
        / records.len() as f64;

    let elapsed = start.elapsed();
    println!("Time elapsed in sequential() is: {:?}", elapsed);
    println!("Average quality score: {}", average);
    average
}

fn mpsc_calculate_average(records: Vec<Record>) -> f64 {
    let start = Instant::now();
    let (tx, rx) = mpsc::channel();

    let record_handle = thread::spawn(move || {
        println!("Hello from the mpsc thread!");
        for record in records {
            let val = calculate_fastq_quality_score(&record.qual);
            tx.send(val).unwrap();
        }
    });

    record_handle.join().unwrap();

    let elapsed = start.elapsed();
    println!("Time elapsed in mpsc_calculate_average() is: {:?}", elapsed);

    return rx.iter().sum::<f64>();
}

fn rayon_calculate_average(records: &Vec<Record>) -> f64 {
    let start = Instant::now();
    let average = records
        .par_iter()
        .map(|r| calculate_fastq_quality_score(&r.qual))
        .sum::<f64>()
        / records.len() as f64;

    let elapsed = start.elapsed();
    println!(
        "Time elapsed in rayon_calculate_average() is: {:?}",
        elapsed
    );
    println!("Average quality score: {}", average);
    average
}

fn main() {
    let fastqfiles = fefastq::read_and_get_records().unwrap();
    for fq in fastqfiles {
        println!("{}", fq.fp);
        one_to_one_mpsc(fq.records.clone());
        sequential(&fq.records);
        rayon_calculate_average(&fq.records);
        let val = mpsc_calculate_average(fq.records);
        println!("Avg Score for {}: {}", fq.fp, val);
    }
}
