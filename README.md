# RNAseq

by Judith Bergadà Pijuan

This RNAseq pipeline is a tool to perform end-to-end RNAseq analyses.
Given sequencing reads (FASTQ files), it provides a list of differentially expressed genes across conditions.

## Installation

To use this pipeline, you need to install the following dependencies:
- STAR >= 2.7.3a
- FastQC >= v0.11.9
- R >= 3.6.0

Later, you need to download the tool:
```bash
cd $HOME
git clone https://github.com/judithbergada/Pipeline_RNAseq
```

## Usage

The pipeline expects you to have the following folders:
- FASTQ folder: this is a folder containg only your sequencing reads (FASTQ files).
- Annotation folder: this is a folder containing only 2 files (a FASTA file and a GTF file).

In addition, you will need at least one text file showing the experimental design.
It should look as follows:
```bash
Rep1_3_S24_R1_001	condition
Rep2_3_S31_R1_001	condition
Rep3_3_S38_R1_001	condition
Rep1_4_S25_R1_001	control
Rep2_4_S32_R1_001	control
Rep3_4_S39_R1_001	control
```

To get information about the usage, please try:
```bash
./RNAseq.sh -h
```

The RNAseq tool can be used with these parameters:
```bash
Usage: RNAseq     [-h or --help]
                  [-f or --fastqfolder]
                  [-a or --annotfolder]
                  [-c or --condfiles]
                  [-o or --outname]
                  [-t or --threads]

Optional arguments:
    -h, --help:
                Show this help message and exit.
    -o, --outname:
                Name of your analysis.
                It will be used to name the output files.
                Default: rnaseq.
    -t, --threads:
                Number of threads that will be used.
                It must be an integer.
                Default: 8.

Required arguments:
    -f, --fastqfolder:
                Path to the folder that contains ALL your FASTQ files.
                Only FASTQ files should be placed in it.
                Files might (or might not) be compressed.
    -a, --annotfolder:
                Path to the folder that contains the annotation files.
                Only 2 files should be placed in this folder:
                1) Whole genome FASTA file, 2) GTF file.
    -c, --condfiles:
                Path to a tab-delimited text file describing the
                experimental design or comparisons. 2 columns are needed:
                Column 1: name of each FASTQ file, without extension.
                Column 2: condition of each file (experiment or control).
                You can specify as many files as you want, but
                they must be writen within quotes.
                E.g.: --condfiles 'comparison1.txt comparison2.txt'
```

Enjoy using the tool!

