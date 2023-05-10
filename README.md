# FeQC: Utilities for FASTQ in Rust
This is a **work in progress**, but there are some useful features implemented.

## FeQC

### Usage

    feqc -i [input/dir/with/fastqs]

The output is a newline for each file and the corresponding average FASTQ quality. For example:

    Staph_R1_head.fastq     33.87047619047619
    Staph_R2_head.fastq     34.24922088040515