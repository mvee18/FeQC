# FeQC: Utilities for FASTQ in Rust
This is a **work in progress**, but there are some useful features implemented.

## FeQC

### Usage

    feqc -i [input/dir/with/fastqs]

The output is a newline for each file and the corresponding average FASTQ quality. For example:

    Staph_R1_head.fastq     33.87047619047619
    Staph_R2_head.fastq     34.24922088040515

Additionally, a new file will be created with the same basename as the fastq files, but with a .qual extensions, that contains a TSV of index and average quality at that position. That is, if `ecoli_R1.fastq` is a read, then a file called `ecoli_R1.qual` will be created in the same directory as the original file. It has the following format:

    Index	Average Quality Score
    0	33
    1	33
    2	32.5