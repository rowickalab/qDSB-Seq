#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

if( @ARGV<3 ){
    print STDERR "perl $0 <btt file> <chr length file: chr length> <output prefix> <ignore zero depth: 0 (default) or 1> <first or full read depth or fragment depth: 1 or 2 or 3>\n

This version will cost more memory but speed improvement

";
    exit;
}

my $btt_file=shift;
my $chrLength=shift;
my $prefix=shift;
my $ignore_zero=shift || 0;
my $first_depth=shift || 1;

my $href_depth={};
if( $first_depth == 1 ){
    $href_depth=btt_to_depth_first_base($btt_file);
}elsif( $first_depth == 2 ){
    $href_depth=btt_to_depth_read($btt_file);
}elsif( $first_depth == 3 ){
    $href_depth=btt_to_depth_fragment($btt_file);
}

my $href_chr=getChrLength($chrLength);

my $chr_dir="";
my ($out1,$out2,$out3,$out4)="";
if( $first_depth == 1 ){
    $chr_dir="$prefix\.startDepth.byChr";
    $out1="$prefix.startDepth";
    $out2="$prefix.startDepth.plus_and_minus";
    $out3="$prefix.startDepth.plus_and_minus_sameline";
    $out4="$prefix.startDepth.plus-minus";
}elsif( $first_depth == 2 ){
    $chr_dir="$prefix\.readDepth.byChr";
    $out1="$prefix.readDepth";
    $out2="$prefix.readDepth.plus_and_minus";
    $out3="$prefix.readDepth.plus_and_minus_sameline";
    $out4="$prefix.readDepth.plus-minus";
}elsif( $first_depth == 3 ){
    $chr_dir="$prefix\.fragDepth.byChr";
    $out1="$prefix.fragDepth";
    $out2="$prefix.fragDepth.plus_and_minus";
    $out3="$prefix.fragDepth.plus_and_minus_sameline";
    $out4="$prefix.fragDepth.plus-minus";
}

`rm -rf $chr_dir $out1.bedGraph $out2.bedGraph $out3.bedGraph $out4.bedGraph`;
mkdir $chr_dir;


foreach my $chr(sort keys %$href_depth){
    print STDERR $chr,"\n";
    my $h=$href_depth->{$chr};

    open(OUT1," > $chr_dir/$out1.$chr.bedGraph");
    open(OUT2," > $chr_dir/$out2.$chr.bedGraph");
    open(OUT4," > $chr_dir/$out3.$chr.bedGraph");
    open(OUT3," > $chr_dir/$out4.$chr.bedGraph");
    my $o1="";
    my $o2="";
    my $o3="";
    my $o4="";

    for( my $p=1;$p<=$href_chr->{$chr};$p++){
        $h->{$p}->{'+'}=0 if !$h->{$p}->{'+'};
        $h->{$p}->{'-'}=0 if !$h->{$p}->{'-'};
        if( $ignore_zero == 0 ){
            $o1.="$chr\t$p\t".($p+1)."\t".($h->{$p}->{'+'}+$h->{$p}->{'-'})."\n";
        }elsif( $ignore_zero == 1 && $h->{$p}->{'+'}+$h->{$p}->{'-'} != 0 ){
            $o1.="$chr\t$p\t".($p+1)."\t".($h->{$p}->{'+'}+$h->{$p}->{'-'})."\n";
        }
       
        if( $ignore_zero == 0 ){
            $o2.="$chr\t$p\t".($p+1)."\t$h->{$p}->{'+'}\n";
            $o2.="$chr\t$p\t".($p+1)."\t-$h->{$p}->{'-'}\n";
            $o4.="$chr\t$p\t".($p+1)."\t$h->{$p}->{'+'}\t-$h->{$p}->{'-'}\n";
        }elsif( $ignore_zero == 1 ){
            if( $h->{$p}->{'+'} != 0 ){
                $o2.="$chr\t$p\t".($p+1)."\t$h->{$p}->{'+'}\n";
                $h->{$p}->{'-'} = 0 if not defined $h->{$p}->{'-'};
            }
            if( $h->{$p}->{'-'} != 0 ){
                $o2.="$chr\t$p\t".($p+1)."\t-$h->{$p}->{'-'}\n";
                $h->{$p}->{'+'} = 0 if not defined $h->{$p}->{'+'};
            }
            if( $h->{$p}->{'+'} != 0 || $h->{$p}->{'-'} != 0 ){
                $o4.="$chr\t$p\t".($p+1)."\t$h->{$p}->{'+'}\t-$h->{$p}->{'-'}\n";
            }
        }

        if( $ignore_zero == 0 ){
            $o3.="$chr\t$p\t".($p+1)."\t".($h->{$p}->{'+'}-$h->{$p}->{'-'})."\n";
        }elsif( $ignore_zero == 1 && $h->{$p}->{'+'}-$h->{$p}->{'-'} != 0 ){
            $o3.="$chr\t$p\t".($p+1)."\t".($h->{$p}->{'+'}-$h->{$p}->{'-'})."\n";
        }
        delete $h->{$p};
    }

    print OUT1 $o1;
    print OUT2 $o2;
    print OUT3 $o3;
    print OUT4 $o4;

   `cat $chr_dir/$out1.$chr.bedGraph >> $out1.bedGraph`;
   `cat $chr_dir/$out2.$chr.bedGraph >> $out2.bedGraph`;
   `cat $chr_dir/$out3.$chr.bedGraph >> $out3.bedGraph`;
   `cat $chr_dir/$out4.$chr.bedGraph >> $out4.bedGraph`;

    close OUT1;
    close OUT2;
    close OUT3;
    close OUT4;
}


sub btt_to_depth_first_base{
    my $file=shift;
    my $href_depth={};

    open(IN,$file);
    while(<IN>){
        chomp;
        next if !$_;
        my @cols=split /\t/;
        if( $cols[1] eq '+' ){ 
            $href_depth->{$cols[2]}->{$cols[3]+1}->{'+'}+=1;
        }elsif( $cols[1] eq '-' ){
            $href_depth->{$cols[2]}->{$cols[3]+1}->{'-'}+=1;
        }
    } 
    close IN;
    return $href_depth;
}

sub btt_to_depth_read{
    my $file=shift;
    my $href_depth={};

    my $n=0;
    open(IN,$file);
    while(<IN>){
        chomp;
        next if !$_;
        my @cols=split /\t/;
        $n++;
        if( $n % 100000 == 0 ){
            print STDERR "$n\n";
        }
        my $seq=$cols[4];
        my $N=length $seq;
        for( my $i=0;$i<$N;$i++ ){
            if( $cols[1] eq '+' ){ 
                my $pos=$cols[3]+1+$i;
                $href_depth->{$cols[2]}->{$pos}->{'+'}+=1;
            }elsif( $cols[1] eq '-' ){
                my $pos=$cols[3]+1-$i;
                $href_depth->{$cols[2]}->{$pos}->{'-'}+=1;
            }
        }
    } 
    close IN;
    return $href_depth;
}

sub btt_to_depth_fragment{
    my $file=shift;
    my $href_depth={};

    my $n=0;
    open(IN,$file);
    while(<IN>){
        chomp;
        next if !$_;
        my @cols_1=split /\t/;

        my $line_2=<IN>;
        my @cols_2=split /\t/,$line_2;

        $n++;
        if( $n % 100000 == 0 ){
            print STDERR "$n\n";
        }

        my ($id_1,$idx_1)=split /\//,$cols_1[0];
        my ($id_2,$idx_2)=split /\//,$cols_2[0];

        if( $cols_1[2] eq $cols_2[2] && $id_1 eq $id_2 ){
            my $start=0;
            my $end=0;
            my $strand='';
            if( $cols_1[1] eq '+' && $cols_2[1] eq '-' ){
                $start=$cols_1[3];
                $end  =$cols_2[3];
                # if R1 is '+', it is '+' strand fragment; if R1 is '-', it is '-' strand fragment
                if( $idx_1 == 1 && $idx_2 == 2 ){
                     $strand='+';
                }elsif( $idx_1 == 2 && $idx_2 == 1 ){
                     $strand='-';
                }else{
                     print STDERR "@cols_1\n";
                     print STDERR "@cols_2\n";
                     exit;
                }
            }elsif( $cols_1[1] eq '-' && $cols_2[1] eq '+' ){
                $start=$cols_2[3];
                $end  =$cols_1[3];
                if( $idx_1 == 1 && $idx_2 == 2 ){
                     $strand='-';
                }elsif( $idx_1 == 2 && $idx_2 == 1 ){
                     $strand='+';
                }else{
                     print STDERR "@cols_1\n";
                     print STDERR "@cols_2\n";
                     exit;
                }
            }else{
                print STDERR "@cols_1\n";
                print STDERR "@cols_2\n";
                exit;
            }

            for( my $pos=$start;$pos<=$end;$pos++ ){
                $href_depth->{$cols_1[2]}->{$pos}->{'+'}+=1;
            }

        }else{
                print STDERR "@cols_1\n";
                print STDERR "@cols_2\n";
                exit;
        }

    } 
    close IN;
    return $href_depth;
}

sub getChrLength{
    my $file=shift;
    my $href={};
    open(IN,$file);
    while(<IN>){
        chomp;
        my ($chr,$length)=split;
        $href->{$chr}=$length;
    }
    close IN;
    return $href;
}
