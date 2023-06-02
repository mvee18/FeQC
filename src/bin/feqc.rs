use fefastq::{average_quality_at_index, rayon_get_average_quality_score};
use std::fs::File;
use std::io::Write;
use std::path::Path;

fn main() {
    let fastqfiles = fefastq::read_and_get_records().unwrap();
    for fq in fastqfiles {
        let avg = rayon_get_average_quality_score(&fq.records);
        let fp = Path::new(&fq.fp);

        println!("{}\t{}", fp.file_name().unwrap().to_str().unwrap(), avg);

        let index_qual = average_quality_at_index(&fq.records);

        let mut file = File::create(fp.with_extension("qual")).unwrap();

        let header = "Index\tAverage Quality Score\n";
        file.write(header.as_bytes()).unwrap();
        for (i, q) in index_qual.iter().enumerate() {
            writeln!(file, "{}\t{}", i, q).unwrap();
        }
    }
}
