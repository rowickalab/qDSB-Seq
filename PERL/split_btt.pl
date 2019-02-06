#!/usr/bin/perl
use strict;
use warnings;
use autodie;

my $btt=shift;

if( ! -d "$btt.split" ){
  mkdir "$btt.split";
}

open my $in_fh, '<', $btt;

my %out_fh;

while (<$in_fh>) {
  next if $_ !~ /\S+/;
  my @cols=split; 
  my $filename = "$btt.split/$cols[2]";
  open $out_fh{$filename}, '>', $filename unless $out_fh{$filename};
  print { $out_fh{$filename} } $_;
}

close $_ for values %out_fh;

