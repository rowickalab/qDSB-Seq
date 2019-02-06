#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

if( @ARGV<2 ){
    print "
This script counts the number of reads in a specific region (eg. G-quadruplex) included in a bed file (3 columns: chr start end).
The sequencing depth for each position on the genome is saved in a wig file (4 colums: chr start end depth). 

Usage: perl $0 <bed_file> <wig_file> [options] 

    -d CHARACTER	output directory
    -p CHARACTER	output prefix
    -l NUMERIC		extend from left position to the left
    -r NUMERIC		extend from right position to the right
    -center		if defined: take center of two postions in bed file, 
    			if not: take two positions in bed file
    -checkovl		if defined: check overlap, and redefine region from the middle of overlapped regions
                        if not: don't check overlap, use left and right poistions 
Output:
1) .txt	record the number of reads in each region
2) .stat	a brief statistics
3) .pos	record the details of reads in each region and position
4) .mat	record the details of reads in each region and position by matrix

";
    exit;
}

my $bed_file=shift; # specific reg where you are intrst
my $depth_file=shift; #plus and minus depth

my $outdir="outdir";
my $outPrefix="output";
my $take_center=0; # take center or not, yes or no
my $flanking_size_left=0;
my $flanking_size_right=0;
my $check_overlap=0;

GetOptions( 'd=s' => \$outdir,
            'p=s' => \$outPrefix,
            'l=i'  => \$flanking_size_left,
            'r=i'  => \$flanking_size_right,
            'checkovl' => \$check_overlap,
            'center' => \$take_center );

my ($href_depth,$n_of_total_mapping_reads,$n_of_total_mapping_reads_plus,$n_of_total_mapping_reads_minus,$href_depth_chr)=loadPlusMinusFile($depth_file);
my $aref_bed=loadBED($bed_file);

mkdir $outdir if !-d $outdir;

open(O1,"> $outdir/$outPrefix.flank_$flanking_size_left\_$flanking_size_right.txt");
open(O2,"> $outdir/$outPrefix.flank_$flanking_size_left\_$flanking_size_right.stat");
open(O3,"> $outdir/$outPrefix.flank_$flanking_size_left\_$flanking_size_right.pos");
open(O4,"> $outdir/$outPrefix.flank_$flanking_size_left\_$flanking_size_right.mat");

# header
print O1 "Chr\tBED_start\tBED_end\tName\tStart\tEnd\tRegion_size\tTotal_reads\tPerc_total\tPlus_reads\tPerc_plus\tMinus_reads\tPerc_minus\tPeak_reads\tPerc_peak_reads\tPlus_peak_reads\tPerc_plus_peak_reads\tMinus_peak_reads\tPerc_minus_peak_reads\tSummit_reads\tPerc_summit_reads\n";

my $n_of_intrst_reg=0;
my $reads_total=0;
my $reads_plus=0;
my $reads_minus=0;
my $reads_peak_depth=0;
my $reads_plus_peak_depth=0;
my $reads_minus_peak_depth=0;
my $reads_summit_depth=0;
my $feat_size=0;

my $n=0;

open(IN,$bed_file);
while(<IN>){
    next if /^#/;
    $n++;
    chomp;
    my ($chr,$s,$e,$name)=split;

    # assign name
    if( !$name ){
        $name="ID$n";
    }

    # compute start and end
    my $center=int(($s+$e)/2);
    my $start=0;
    my $end  =0;
    if( $take_center == 1 ){
        $start=$center-$flanking_size_left+1;
        $end  =$center+$flanking_size_right;
    }elsif( $take_center == 0 ){
        $start=$s-$flanking_size_left;
        $end  =$e+$flanking_size_right;
    }else{
        print STDERR "Please check your 'center' option\n";
        exit;
    }

    if( $start < 1 ){
        $start=1;
    }

    if( $check_overlap == 1 ){
      # check with previous
      if( $aref_bed->[$n-2] ){
        my $aref_previous=$aref_bed->[$n-2];
        my $previous_chr=$aref_previous->[0];
        my $previous_center=int(($aref_previous->[1]+$aref_previous->[2])/2);
        if( $chr eq $previous_chr ){
            my $previous_end=$previous_center+$flanking_size_right;
            if( $start - $previous_end < 0 ){
                $start = int(($center+$previous_center)/2);
            }
        }
      }

      # check with next
      if( $aref_bed->[$n] ){
        my $aref_next=$aref_bed->[$n];
        my $next_chr=$aref_next->[0];
        my $next_center=int(($aref_next->[1]+$aref_next->[2])/2);
        if( $chr eq $next_chr ){
            my $next_start=$next_center-$flanking_size_left;
            if( $end - $next_start > 0 ){
                $end = int(($center+$next_center)/2);
            }
        }
      }
    }

    my ($plus,$minus,$total,$region_size,$density,$plus_peak_depth,$minus_peak_depth,$depth_for_each_position,$depth)=getDepth($href_depth,$chr,$start,$end);

    my $peak_depth=$plus_peak_depth+$minus_peak_depth;
    
    my $summit_depth=0;
    if( $plus_peak_depth > $minus_peak_depth ){
        $summit_depth=$plus_peak_depth;
    }else{
        $summit_depth=$minus_peak_depth;
    }

    my $perc_total=sprintf("%.2f",$total/$n_of_total_mapping_reads*100);
    my $perc_plus=sprintf("%.2f",$plus/$n_of_total_mapping_reads_plus*100);
    my $perc_minus=sprintf("%.2f",$minus/$n_of_total_mapping_reads_minus*100);

    my $perc_peak_depth=sprintf("%.2f",$peak_depth/$n_of_total_mapping_reads*100);
    my $perc_summit_depth=sprintf("%.2f",$summit_depth/$n_of_total_mapping_reads*100);
    my $perc_plus_peak_depth=sprintf("%.2f",$plus_peak_depth/$n_of_total_mapping_reads_plus*100);
    my $perc_minus_peak_depth=sprintf("%.2f",$minus_peak_depth/$n_of_total_mapping_reads_minus*100);

    print O1 "$chr\t$s\t$e\t$name\t$start\t$end\t$region_size\t$total\t$perc_total\t$plus\t$perc_plus\t$minus\t$perc_minus\t$peak_depth\t$perc_peak_depth\t$plus_peak_depth\t$perc_plus_peak_depth\t$minus_peak_depth\t$perc_minus_peak_depth\t$summit_depth\t$perc_summit_depth\n";
    print O3 "#$chr\t$start\t$end\n$depth_for_each_position";
    print O4 "$chr\t$start\t$end\t",join("\t",@$depth),"\n";

    $n_of_intrst_reg++;
    $reads_total+=$total;
    $reads_plus+=$plus;
    $reads_minus+=$minus;
    $reads_peak_depth+=$peak_depth;
    $reads_plus_peak_depth+=$plus_peak_depth;
    $reads_minus_peak_depth+=$minus_peak_depth;
    $reads_summit_depth+=$summit_depth;
    $feat_size+=$end-$start+1;
}

my $perc_reads_total=sprintf("%.2f",$reads_total/$n_of_total_mapping_reads*100);
my $perc_reads_plus=sprintf("%.2f",$reads_plus/$n_of_total_mapping_reads_plus*100);
my $perc_reads_minus=sprintf("%.2f",$reads_minus/$n_of_total_mapping_reads_minus*100);
my $perc_peak_depth=sprintf("%.2f",$reads_peak_depth/$n_of_total_mapping_reads*100);
my $perc_plus_peak_depth=sprintf("%.2f",$reads_plus_peak_depth/$n_of_total_mapping_reads_plus*100);
my $perc_minus_peak_depth=sprintf("%.2f",$reads_minus_peak_depth/$n_of_total_mapping_reads_minus*100);
my $perc_summit_depth=sprintf("%.2f",$reads_summit_depth/$n_of_total_mapping_reads*100);

print O2 "Summary
Total reads:\t$n_of_total_mapping_reads
Plus  reads:\t$n_of_total_mapping_reads_plus
Minus reads:\t$n_of_total_mapping_reads_minus
Features:\t$n_of_intrst_reg
Feature size:\t$feat_size
Feature total reads:\t$reads_total\t$perc_reads_total%
Feature plus  reads:\t$reads_plus\t$perc_reads_plus%
Feature minus reads:\t$reads_minus\t$perc_reads_minus%
Feature peak  reads:\t$reads_peak_depth\t$perc_peak_depth%
Feature +peak reads:\t$reads_plus_peak_depth\t$perc_plus_peak_depth%
Feature -peak reads:\t$reads_minus_peak_depth\t$perc_minus_peak_depth%
Feature summit reads:\t$reads_summit_depth\t$perc_summit_depth%

";


sub loadPlusMinusFile{
    my $file=shift;
    
    my $href_depth={};
    my $total_depth=0;
    my $total_depth_plus=0;
    my $total_depth_minus=0;
    my $href_depth_chr={};
    
    open(IN,$file);
    while(<IN>){
        chomp;
        my @cols=split;
        if( $cols[3] > 0 ){
            $href_depth->{$cols[0]}->{$cols[1]}->{'plus'}=abs($cols[3]);
            $total_depth_plus+=abs($cols[3]);
            $href_depth_chr->{$cols[0]}+=abs($cols[3]);
        }elsif( $cols[3] < 0 ){
            $href_depth->{$cols[0]}->{$cols[1]}->{'minus'}=abs($cols[3]);
            $total_depth_minus+=abs($cols[3]);
            $href_depth_chr->{$cols[0]}+=abs($cols[3]);
        }
    }
    close IN;

    $total_depth=$total_depth_plus+$total_depth_minus;
    return $href_depth,$total_depth,$total_depth_plus,$total_depth_minus,$href_depth_chr;
}

sub getDepth{
    my ($href_depth,$chr,$start,$end)=@_;

    my ($plus,$minus,$plus_peak_depth,$minus_peak_depth,$depth_for_each_position)=0;
    my @depth;
    if( not defined $plus_peak_depth ){
        $plus_peak_depth=0;
    }
    if( not defined $minus_peak_depth ){
        $minus_peak_depth=0;
    }

    my $relative_pos_to_center=0;
    my $center=int(abs($start+$end)/2);
    for( my $pos=$start;$pos<=$end;$pos++){
        my ($p,$m,$t)=0;
        my $rel=$pos-$center;
        if( defined $href_depth->{$chr}->{$pos}->{'plus'} ){
            $p=$href_depth->{$chr}->{$pos}->{'plus'};
            $plus+=$href_depth->{$chr}->{$pos}->{'plus'};
            if( $href_depth->{$chr}->{$pos}->{'plus'} > $plus_peak_depth ){
                $plus_peak_depth=$href_depth->{$chr}->{$pos}->{'plus'};
            }
        }
        if( defined $href_depth->{$chr}->{$pos}->{'minus'} ){
            $m=$href_depth->{$chr}->{$pos}->{'minus'};
            $minus+=$href_depth->{$chr}->{$pos}->{'minus'};
            if( $href_depth->{$chr}->{$pos}->{'minus'} > $minus_peak_depth ){
                $minus_peak_depth=$href_depth->{$chr}->{$pos}->{'minus'};
            }
        }
     
        # 
        $p = 0 if !$p;
        $m = 0 if !$m;
        $t = $p + $m;
        $depth_for_each_position.="$chr\t$pos\t$rel\t$t\t$p\t$m\n";
        push @depth, $t;
    }
    $plus = 0 if !$plus;
    $minus= 0 if !$minus;
    $plus_peak_depth=0 if !$plus_peak_depth;
    $minus_peak_depth=0 if !$minus_peak_depth;

    my $total=abs($plus)+abs($minus);
    my $region_size=$end-$start+1;
    print "$chr,$start,$end\n";
    my $density=$total/$region_size;

    return abs($plus),abs($minus),$total,$region_size,$density,abs($plus_peak_depth),abs($minus_peak_depth),$depth_for_each_position,\@depth;
}

sub loadBED{
    my $file=shift;
    my @array;
    open(IN,$file);
    while(<IN>){
        next if /^#/;
        chomp;
        my @cols=split;
        push @array, \@cols;
    }
    return \@array;
}

close O1;
close O2;
close O3;
close O4;
