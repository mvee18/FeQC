use fefastq::rayon_get_average_quality_score;
use std::path::Path;

fn main() {
    let fastqfiles = fefastq::read_and_get_records().unwrap();
    for fq in fastqfiles {
        let avg = rayon_get_average_quality_score(&fq.records);
        let fp = Path::new(&fq.fp);

        println!("{}\t{}", fp.file_name().unwrap().to_str().unwrap(), avg);
    }
}
