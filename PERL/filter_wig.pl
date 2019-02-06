#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

if( @ARGV<4 ){
    print STDERR "Usage: perl $0 <wig file> <flanking size of bed> <ignore locations in bed: 1 or not 0> <bed files>\n";
    exit;
}

my $wig_file=shift;
my $flanking_size=shift;
my $ignore=shift;
my @bed_files=@ARGV;

my $href_bed={};
foreach my $bed_file(@bed_files){
    open(BED,$bed_file);
    while(<BED>){
        chomp;
        my @cols=split;
        push @{$href_bed->{$cols[0]}}, [$cols[0],$cols[1]-$flanking_size,$cols[2]+$flanking_size];
    }
    close BED;
}

open(WIG,$wig_file);
while(<WIG>){

    my @cols=split;

    next if !$cols[0];
    next if $cols[3] == 0;

    my $bed=$href_bed->{$cols[0]};

    if( $ignore == 1 ){
        if( search($cols[0],$cols[1],$bed) == 0 ){
            print;
        }
    }else{
        if( search($cols[0],$cols[1],$bed) == 1 ){
            print;
        }
    }
}

sub search{
    my $q_chr=shift;
    my $q_start=shift;
    my $aref=shift;

    foreach my $a(@$aref){
        if( $q_chr eq $a->[0] ){
            if( $q_start >= $a->[1] && $q_start <= $a->[2] ){
                return 1;
            }
        }
    }

    return 0;
}
