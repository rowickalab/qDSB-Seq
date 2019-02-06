#!/usr/bin/perl -w
use strict;
use warnings;

if(@ARGV<1){
	print "Usage: perl $0 <fasta>\n";
	exit;
}

my $input=shift;

if( -d $input ){
	opendir(DIR,$input);
	my @files=grep /(\.fa.?.?.?|gb)/,readdir DIR;
	foreach my $file(@files){
		fasLength("$input/$file");
	}
}else{
	fasLength("$input");
}

######################################
sub fasLength{
    my $fasta = shift;
    my $name = '';
    my $seq = '';
    my $length = 0;
    open(IN,"$fasta") or die "Can't open $!";
    while (<IN>){
        chomp;
	s/\r$//;
        if(/^>(\S+)/){
	    if( $seq ){
                $length += length($seq);
	        $name=~s/\s.*//;
	        print "$name\t$length\n";
	        print STDERR "$name\t1\t$length\n";
                $seq = ();
	    }
            $name = $1;
	    $length = 0;
        }elsif( eof ){
	    s/\s//g;
            $seq .= $_;
	    if( $seq ){
                $length += length($seq);
	        $name=~s/\s.*//;
	        print "$name\t$length\n";
	        print STDERR "$name\t1\t$length\n";
                $seq = ();
	    }
        }elsif(/^[A-Za-z]/){
	    s/\s//g;
            $seq .= $_;
        }
    }
    #print "base pairs of all sequences: $length_sum\n";
    close IN;
}

