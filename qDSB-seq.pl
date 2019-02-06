#!/usr/bin/perl
use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Long;

if( @ARGV<2 ){
    print "
qDSB-seq.pl performs quantitative DSB sequencing analysis, including:

1) map DSB sequencing data to genome, and call depth
2) map gDNA sequencing data to genome, and calculate enzyme cutting efficiency
3) calculate DSB frequencies per cell on the whole genome or specific locations 

Usage: perl $0 <DSB reads> <gDNA reads R1> <gDNA reads R2> [options] 

    -s CHARACTER	sample name
    -r CHARACTER        group of the sample, for example, G1 or S
    -g CHARACTER	genome name
    -f CHARACTER        genome sequence in fasta format
    -i CHARACTER	genome bowtie index
    -e CHARACTER	enzyme name
    -t CHARACTER        enzyme type
    -c CHARACTER        enzyme cutting sites in BED format
    -b CHARACTER        backgrond coordinates on genome to remove background noise from cutting efficiency calculation
    -p CHARACTER        output prefix

Please contact Dr. Yingjie Zhu for any questions on the code yizhu\@utmb.edu

";
    exit;
}

my $DSB_reads=shift;
my $gDNA_reads_1=shift;
my $gDNA_reads_2=shift;

my $sample_name="unknown";
my $group_sample="unknown";
my $genome_name="genome";
my $bowtie_index="";
my $genome_sequence="";
my $enzyme_name="unknown";
my $enzyme_type="3";
my $enzyme_cutting_sites="";
my $background="";
my $prefix="output";

GetOptions( 's=s' => \$sample_name,
            'r=s' => \$group_sample,
            'g=s' => \$genome_name,
            'f=s' => \$genome_sequence,
            'i=s' => \$bowtie_index,
            'e=s' => \$enzyme_name,
            't=s' => \$enzyme_type,
            'c=s' => \$enzyme_cutting_sites,
            'b=s' => \$background,
            'p=s' => \$prefix
             );

my $chr_length="$genome_sequence.length";
if( $genome_sequence ){
    `$Bin/PERL/fastaLength.pl $genome_sequence > $genome_sequence.length 2>$genome_sequence.length.bed`;
}

`sh $Bin/process_DSB-seq_data.sh $genome_name $chr_length $bowtie_index $enzyme_name $enzyme_type $enzyme_cutting_sites $DSB_reads $prefix`;

`sh $Bin/process_gDNA_data.sh $genome_name $chr_length $bowtie_index $enzyme_name $enzyme_type $enzyme_cutting_sites $gDNA_reads_1 $gDNA_reads_2 $background $prefix`;

my $switch='';
if( $enzyme_type eq '5' ){
    $switch='.switch';
}
    
`Rscript $Bin/qDSB-seq.R -s $sample_name -e $enzyme_name -g $group_sample -m process_gDNA_data/cutting_efficiency_$enzyme_name/$prefix.$enzyme_name.L0-R0.exact$switch.txt -b process_gDNA_data/cutting_efficiency_$enzyme_name/$prefix.$enzyme_name.background.exact$switch.txt -r process_DSB-seq_data/get_$enzyme_name\_reads/$prefix.$enzyme_name.flank_0_0.txt -d process_DSB-seq_data/map_SE_to_$genome_name\_btt_depth/$prefix.startDepth.plus_and_minus.wig -t 2 -c 10 -n 0 -a 1 -p $prefix`;
