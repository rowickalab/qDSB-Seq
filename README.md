# qDSB-Seq
Custom code for quantitative DSB sequencing (qDSB-Seq)

# Overview
With the development of label-based DSB sequencing technologies, DSB quantification and normalization is becoming important. qDSB-Seq, a quantitative DSB sequencing strategy, provides an optimized solution to measure DSB numbers per cell and their precise genomic coordinates at the same time. Here we provide a package for qDSB-Seq analysis. 

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

    git clone https://github.com/rowickalab/qDSB-Seq.git

Or download the package and then unzip:

    unzip qDSB-Seq-master.zip

Compile btt software that convert bowtie output (gcc required):

    cd src
    make

Bowtie should be installed before running the code:
1) download bowtie: http://bowtie-bio.sourceforge.net/index.shtml
2) install it to the local directory
3) add bowtie software to the path that can be found by enviroment

Please download hygestat_bless in iSeq package for read preprocessing: https://github.com/rowickalab/iSeq

# Instruction for use

qDSB-Seq.pl integrates the scripts written by R, PERL, C++, and BASH for an easy use. To get help, type ‘perl qDSB-Seq.pl’ on the command line of Linux. 

qDSB-seq.pl performs quantitative DSB sequencing analysis, including:

    1) map DSB sequencing data to genome, and call depth
    2) map gDNA sequencing data to genome, and calculate enzyme cutting efficiency
    3) calculate DSB frequencies per cell on the whole genome or specific locations

    Usage: perl ../qDSB-seq.pl <DSB reads> <gDNA reads R1> <gDNA reads R2> [options]

    -s CHARACTER        sample name
    -r CHARACTER        group of the sample, for example, G1 or S
    -g CHARACTER        genome name
    -f CHARACTER        genome sequence in fasta format
    -i CHARACTER        genome bowtie index
    -e CHARACTER        enzyme name
    -t CHARACTER        enzyme type
    -c CHARACTER        enzyme cutting sites in BED format
    -b CHARACTER        backgrond coordinates on genome to remove background noise from cutting efficiency calculation
    -p CHARACTER        output prefix

Here who can follow the example to learn how to run this code. This example come from real data, but a selected dataset. The genome is cleaved by NotI enzyme.

Before running the code, who should prepare or use the input data as follows:
 
1.Sequencing reads from DSB sequencing, only sequence without name and quality inside
    
    test_i-BLESS.seq
    
    GGCCGCCACCATCGCGATGGTAACGGCAGTAGCAACGGTAATGGTGAAC
    GGCCGCCACCATCGCGATGGTAACGGCAGTAGCAACGGTAATGGTGAACC
    GGCCGCCACCATCGCGATGGTAACGGCAGTAGCAACGGTAATGGTGAAC

2.Paired-end sequencing reads from gDNA sequencing, including R1 and R2 reads, only sequence without name and quality inside

    test_gDNA.R1.seq
    test_gDNA.R2.seq

3.bowtie index of reference genome built by bowtie

    reference_genome/test.reference_genome.bowtie

To build your bowtie index:

    bowtie-build reference_genome index_prefix

4.enzyme cutting sites, it can be obtain from XXX
  
    NotI.bed

5.genome background to remove sequencing fragmentation noise. It can be obtained by running R code on command line
  
    background.bed
   
For your own genome background:

    Rscript ../R/generate_background_bed.R chromosome_length number_of_locations output
    
To run the code for the example:

    cd example

    perl ../qDSB-seq.pl test_i-BLESS.seq test_gDNA.R1.seq test_gDNA.R2.seq -s test -r G1 -g yeast -f reference_genome/test.reference_genome.fas -i reference_genome/test.reference_genome.bowtie -e NotI -t 5 -c NotI.bed -b background.bed -p output
    
The example was cut by NotI enzyme, which cleaves the substrate sequences and creates 5’-overhang. Therefore, we should tell what kind of DNA ends is produced by '-t 5'.

By running this code, it will produce two directories and one summary file. 

    process_DSB-seq_data
    process_gDNA_data
    output.DSBs.summary.txt

process_DSB-seq_data includes mapping file and coverage files on DSB sequencing data. 
process_gDNA_data includes mapping file and cutting efficiency files. Enzyme cutting efficiencies and genome background are under cutting_efficiency_NotI.
output.DSBs.summary.txt contains DSBs per cell on the whole reference genome. 

    sample_name     group   enzyme_name     reads_from_enzyme       reads_from_studied_DSBs fcut_rmbg       DSBs    DSBs_sd_by_site DSBs_sd_by_fcut nsites  cor_fcut_labeled_reads
    test    G1      NotI    498320  501485  0.103111        0.207532        0.06483 0.175323        2       1

# Citation

Please cite our paper on BioRxiv:
Zhu, Y. et al. qDSB-Seq: quantitative DNA double-strand break sequencing. BioRxiv, doi:https://doi.org/10.1101/171405 (2019). 
