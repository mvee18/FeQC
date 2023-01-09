use fefastq::FastqFile;
use std::path::Path;

fn trim_records(fq: FastqFile) {
    let pth = Path::new(&fq.fp);
    let parent = pth.parent().unwrap();
    let stem = pth.file_stem().unwrap();

    // The new name should be parent + stem + _trimmed.fastq
    let new_name = format!("{}_replaced.fastq", stem.to_str().unwrap());

    let new_fp = parent.join(new_name);

    for mut r in fq.records {
        r.plus = String::from("+");
        r.write_record(&new_fp).unwrap();
        println!("{:?}", r)
    }
}
fn main() {
    let fastqfiles = fefastq::read_and_get_records().unwrap();
    for fq in fastqfiles {
        trim_records(fq);
    }
}
