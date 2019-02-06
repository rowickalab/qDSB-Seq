# qDSB-Seq
custom code for quantitative DSB sequencing (qDSB-Seq)

# Overview
With the development of label-based DSB sequencing technologies, DSB quantification and normalization is becoming important. qDSB-Seq, a quantitatively DSB sequencing strategy, provides an optimized solution to measure DSB numbers per cell and their precise genomic coordinates at the same time. Here we provide a package for qDSB-Seq analysis. 

# System requirement
## Hardware Requirements
The qDSB-Seq package requires a computer with enough RAM to support the operations defined by a user. The RAM depends on how big the datasets and the genome size are. For minimal performance, it will be a standard computer with about 8 GB of RAM for our example provided. For optimal performance, we recommend a server with the following specs:

RAM: 256+ GB

CPU: 8+ cores, 3.1+ GHz/core

## Software Requirements
### OS Requirements
This package is supported for Linux operating systems. The package has been tested on the following systems:
Linux: Fedora 20

### Installing R version 3.5.1
1.	Download R from http://cran.us.r-project.org/, click “Download R for Linux” to download latest version.
2.	Install R. Leave all default settings in the installation options.
### R package dependencies 
Once R is installed, type ‘R’ to enter into console, install the packages needed:

    install.packages(“optparse”)

The versions of packages are:
optparse: 1.6.0
### Other requirements 
PERL version 5.18
C++(gcc) compiler

# Installation Guide
To install the package, use git clone:

Or download the package and then unzip:

compile btt software that convert bowtie output (gcc required):

# Instruction for use
qDSB-Seq.pl integrates the scripts written by R, PERL, C++, and BASH for an easy use. To get help, type ‘perl qDSB-Seq.pl’ on the command line of Linux. Here who can follow the example to learn how to use. This example come from real data, but a selected dataset. The genome is cleaved by NotI enzyme. Let’s enter into the ‘example’ directory.
Before running the code, the users should prepare or know the input data as follows:
1)	Sequencing reads from DSB sequencing, only sequence without name and quality inside

    test_i-BLESS.seq
    
    GGCCGCCACCATCGCGATGGTAACGGCAGTAGCAACGGTAATGGTGAAC
    GGCCGCCACCATCGCGATGGTAACGGCAGTAGCAACGGTAATGGTGAACC
    GGCCGCCACCATCGCGATGGTAACGGCAGTAGCAACGGTAATGGTGAAC

2)	Paired-end sequencing reads from gDNA sequencing, including R1 and R2 reads, only sequence without name and quality inside

    test_gDNA.R1.seq
    test_gDNA.R2.seq

3)	bowtie index of reference genome built by bowtie
    reference_genome/test.reference_genome.bowtie

4)	enzyme cutting sites, it can be obtain from XXX
    NotI.bed

5)	genome background to remove sequencing fragmentation noise. It can be obtained by running R code on command line
    background.bed
    
The example was cut by NotI enzyme, which cleaves the substrate sequences and creates 5’-overhang. Therefore, we should tell what kind of DNA ends is produced.

By running this code, it will produce two directories and one summary file. 
For the directory process_DSB-seq_data, it includes mapping file and coverage files on DSB sequencing data. 
For the directory process_gDNA_data, it includes mapping file and cutting efficiency files. Enzyme cutting efficiencies and genome background are under cutting_efficiency_NotI.

# Citation
Please cite our paper on BioRxiv:
Zhu, Y. et al. qDSB-Seq: quantitative DNA double-strand break sequencing. BioRxiv, doi:https://doi.org/10.1101/171405 (2019). 
