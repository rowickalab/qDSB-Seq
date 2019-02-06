#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

if( @ARGV < 2 ){
    print STDERR "
    perl $0 <bed file> <btt file> 

    bed_file    enzyme cutting sites

    btt_file    paired-end mapping output of gDNA sequencing reads from bowtie with uniquely mapping option '-m 1', then output should be converted to btt file using our btt tool

    -rmdup      remove duplicates if it is set

    -l <int>    extend start position to left, default is 0

    -r <int>    extend end position to right, default is 0

    -exact      classify fragment types based on exact position if it is set

    -switch     switch start and end for cutting sites if it is set
                for the cutting site with 5' overhang, left fragment will end with 'end position', right fragment will start with 'start position'

    -summit     compute cutting efficiency based on left or right cut reads with largest value
                cutting efficiency = cut / (cut + uncut reads)

    -prefix     output prefix, default is 'output'

    Please contact Dr. Yingjie Zhu for any questions on the code yizhu\@utmb.edu

\n";
    exit;
}

my $bed_file=shift; # known break coordinates
my $btt_file=shift; # mapping results with btt format

my $remove_duplicate=0; # whether remove duplicates
my $extend_left=0; # how many base pairs extend?
my $extend_right=0;
my $use_exact_position=0; # classify fragment types based on exact position
# compute break efficiency based on summit of LR and RR, 
# break efficiency = summit / (summit + jump).
# Default will use ((LR + RR) / 2) / ((LR + RR) / 2 + jump)
my $use_summit=0; 
my $switch_left_right=0;
my $prefix="output";

GetOptions(
    'rmdup' => \$remove_duplicate,
    'l=i'   => \$extend_left,
    'r=i'   => \$extend_right,
    'exact' => \$use_exact_position,
    'switch' => \$switch_left_right,
    'summit'=> \$use_summit,
    'prefix=s'=> \$prefix
);

# output name
my $output="$prefix.L$extend_left-R$extend_right";
if( $remove_duplicate == 1 ){
    $output.=".rmdup";
}
if( $use_exact_position == 1 ){
    $output.=".exact";
}
if( $use_summit == 1 ){
    $output.=".summit";
}
if( $switch_left_right == 1 ){
    $output.=".switch";
}

# 
open(OUT,"> $output.txt");

print OUT "Chr\tStart\tEnd\tName\tTotal_reads\tUncut_reads\tLeft_reads\tRight_reads\tSummit_reads\tCutting_efficiency\tMin\tMax\n";

# 
my $href_bed=loadBED_byChr($bed_file);


#
foreach my $chr(sort keys %$href_bed){
    print STDERR "$chr\n";
    my $href_frag=btt_fragment_optimized($btt_file,$remove_duplicate,$chr);

    my $aref=$href_bed->{$chr};
    foreach my $aref_bed(@$aref){
        my ($chr,$start,$end,$name)=@$aref_bed;

        if( !$name ){
            $name='-';
        }

        if( $extend_left ){
            $start=$start-$extend_left; }
        if( $extend_right ){
            $end=$end+$extend_right; }

        if( $switch_left_right == 1 ){
            my $tmp=$start;
            $start=$end;
            $end=$tmp;
        } 

        my ($class_results,$href_count_types)=classify_optimized($chr,$start,$end,$href_frag,$use_exact_position);

        my $JR=$href_count_types->{'JR'};
        my $LR=$href_count_types->{'LR'};
        my $RR=$href_count_types->{'RR'};
        my $TR=$JR+$LR+$RR;
        
        my $Summit_reads=0;
        if( $LR > $RR ){
            $Summit_reads=$LR;
        }else{
            $Summit_reads=$RR;
        }

        my ($cutting_frequency,$cutting_frequency_min,$cutting_frequency_max)=calculate_cutting_frequency($use_summit,$TR,$LR,$RR,$JR,$Summit_reads);

        print OUT "$chr\t$start\t$end\t$name\t$TR\t$JR\t$LR\t$RR\t$Summit_reads\t$cutting_frequency\t$cutting_frequency_min\t$cutting_frequency_max\n";
    }
    $href_frag={};
}

close OUT;

sub calculate_cutting_frequency{
    my ($use_summit,$TR,$LR,$RR,$JR,$Summit_reads)=@_;
    #print "$use_summit,$TR,$LR,$RR,$JR,$Summit_reads\n";

    my ($cutting_frequency,$cutting_frequency_min,$cutting_frequency_max)=0;

        if( $use_summit == 1 ){
            my $cut=$Summit_reads;
            my $uncut=$JR;
            if( $cut != 0 ){
                $cutting_frequency = sprintf("%.4f",$cut/($cut+$uncut)) ;
                $cutting_frequency_min = sprintf("%.4f",( $cut-sqrt($cut) ) / ( ($cut-sqrt($cut))+($uncut+sqrt($uncut))) );
                $cutting_frequency_max = sprintf("%.4f",( $cut+sqrt($cut) ) / ( ($cut+sqrt($cut))+($uncut-sqrt($uncut))) );
            }else{
                $cutting_frequency=0;
                $cutting_frequency_min=0;
                $cutting_frequency_max=0;
            }
        }else{
            my $cut=$LR+$RR;
            my $uncut=$JR;
            if( $cut != 0 ){
                $cutting_frequency = sprintf("%.4f",$cut/($cut+2*$JR)) ;
                if( ($cut==1 && $uncut==0 ) ){
                    $cutting_frequency_min=0;
                }else{
                    $cutting_frequency_min = sprintf("%.4f",($cut-sqrt($cut))/(($cut-sqrt($cut))+2*($uncut+sqrt($uncut))) );
                }
                $cutting_frequency_max = sprintf("%.4f",($cut+sqrt($cut))/(($cut+sqrt($cut))+2*($uncut-sqrt($uncut))) );
            }else{
                $cutting_frequency=0;
                $cutting_frequency_min=0;
                $cutting_frequency_max=0;
            }
        } 

    return $cutting_frequency,$cutting_frequency_min,$cutting_frequency_max;

}

sub loadBED_byChr{
 
    my $bed_file=shift;
    my $href={};
    open(IN,"$bed_file");
    while(<IN>){
        chomp;
        next if /^#/;
        my @cols=split /\s+/;
        push @{$href->{$cols[0]}},\@cols;
    }
    return $href;
}

sub classify{
    my $chr=shift;
    my $start=shift;
    my $end=shift;
    my $href_frag=shift;

    my @class_results;
    my $href_count_types={};
    my $href=$href_frag->{$chr};

    foreach my $f_s( sort{$a<=>$b} keys %$href ){
        my $aref_e=$href->{$f_s};
        foreach my $f_e( sort{$a<=>$b} @$aref_e ){
            my $type=classify_fragment($f_s,$f_e,$start,$end); 
            next if $type eq 'NO';
            push @class_results, "$chr\t$f_s\t$f_e\t$type";
            $href_count_types->{$type}++;
        }
    }
    #print join("\n",@class_results),"\n";  
    return \@class_results,$href_count_types;
}

sub classify_optimized{ # using index to improve speed
    my $chr=shift;
    my $start=shift;
    my $end=shift;
    my $href_frag=shift;
    my $use_exact_position=shift; # 1 or 0

    my @class_results;
    my $href_count_types={};
    my $href=$href_frag->{$chr};

    foreach my $idx( keys %$href ){
        my ($idx_start,$idx_end)=split /\t/,$idx;
        if( ( $start >= $idx_start && $start <= $idx_end ) || ( $end >= $idx_start && $end <= $idx_end ) ){
            my $aref_frag=$href->{$idx};
            foreach my $aref( @$aref_frag ){
                my $f_s=$aref->[0];
                my $f_e=$aref->[1];
                my $type='';
                if( $use_exact_position == 1 ){
                    $type=classify_fragment_exact($f_s,$f_e,$start,$end); 
                }else{
                    $type=classify_fragment($f_s,$f_e,$start,$end); 
                }
                next if $type eq 'NO';
                push @class_results, "$chr\t$f_s\t$f_e\t$type";
                $href_count_types->{$type}++;
            }
        }
    }

    if( not defined $href_count_types->{'LR'} ){ $href_count_types->{'LR'} = 0 };
    if( not defined $href_count_types->{'RR'} ){ $href_count_types->{'RR'} = 0 };
    if( not defined $href_count_types->{'JR' } ){ $href_count_types->{'JR' } = 0 };
  
    #print join("\n",@class_results),"\n"; # print read statues  

    return \@class_results,$href_count_types;
}

# classify fragments based on a window 
sub classify_fragment{
    my $f_s=shift;
    my $f_e=shift;
    my $b_s=shift;
    my $b_e=shift;

    if( $f_s >= $f_e ){
        return 'NO'; # the fragment is no relaionship with break
    }elsif( $f_e < $b_s || $f_s > $b_e ){
        return 'NO'; # no overlap between fragment and break
    }elsif( $f_s < $b_s && $f_e > $b_e ){
        return 'JR'; # the fragment cover break
    }elsif( $f_e >= $b_s && $f_e <= $b_e ){
        return 'LR'; # the fragment on the left of break
    }elsif( $f_s >= $b_s && $f_s <= $b_e ){
        return 'RR'; # the fragment on the right of break
    }else{
        return 'Not classified';
    }
}

# classify fragments based on exact position 
sub classify_fragment_exact{
    my $f_s=shift;
    my $f_e=shift;
    my $b_s=shift;
    my $b_e=shift;

    if( $f_s >= $f_e ){
        return 'NO'; # the fragment is no relaionship with break
    }elsif( $f_e < $b_s || $f_s > $b_e ){
        return 'NO';
    }elsif( $f_e == $b_s && $f_s < $b_s ){
        return 'LR'; # the fragment on the left of break
    }elsif( $f_s == $b_e && $f_e > $b_e ){
        return 'RR'; # the fragment on the right of break
    }else{
        return 'JR'; # the fragment cover break
    }
}

sub btt_fragment{
    my $file=shift;
    my $href_frag=shift;

    my $n=0;
    open(IN,$file);
    while(<IN>){
        chomp;
        next if !$_;
        my @cols_1=split /\t/;

        my $line_2=<IN>;
        my @cols_2=split /\t/,$line_2;

        $n++;
        if( $n % 1000000 == 0 ){
            print STDERR "$n\n";
        }

        my ($id_1,$idx_1)=split /\//,$cols_1[0];
        my ($id_2,$idx_2)=split /\//,$cols_2[0];

        if( $cols_1[2] eq $cols_2[2] && $id_1 eq $id_2 ){
            my $start=0;
            my $end=0;
            if( $cols_1[1] eq '+' && $cols_2[1] eq '-' ){
                $start=$cols_1[3];
                $end  =$cols_2[3];
            }elsif( $cols_1[1] eq '-' && $cols_2[1] eq '+' ){
                $start=$cols_2[3];
                $end  =$cols_1[3];
            }else{
                #print STDERR "@cols_1\n";
                #print STDERR "@cols_2\n";
                exit;
            }

            push @{$href_frag->{$cols_1[2]}->{$start}},$end;

        }else{
                #print STDERR "@cols_1\n";
                #print STDERR "@cols_2\n";
                exit;
        }

    } 
    close IN;
    return $href_frag;
}

sub btt_fragment_optimized{ # using index to improve speed
    my $file=shift;
    my $remove_duplicate=shift; # 1 or 0
    my $chr=shift;

    my $href_frag={};
    my $href_check_duplicate={};

    my $n=0;
    open(IN,"$file.split/$chr") || die "Can't find $file.split/$chr, please check you have split directory for btt";
    while(<IN>){
        chomp;
        next if !$_;
        my @cols_1=split /\t/;

        my $line_2=<IN>;
        my @cols_2=split /\t/,$line_2;

        # show progress
        $n++;
        if( $n % 100000 == 0 ){
            print STDERR "$n\n";
        }

        # get read id and read index for R1 and R2
        my ($id_1,$idx_1)=split /\//,$cols_1[0];
        my ($id_2,$idx_2)=split /\//,$cols_2[0];

        # get start and end for paired-end reads
        if( $cols_1[2] eq $cols_2[2] && $id_1 eq $id_2 && $cols_1[2] eq $chr ){
            my $start=0;
            my $end=0;
            if( $cols_1[1] eq '+' && $cols_2[1] eq '-' ){
                $start=$cols_1[3]+1;
                $end  =$cols_2[3]+1;
            }elsif( $cols_1[1] eq '-' && $cols_2[1] eq '+' ){
                $start=$cols_2[3]+1;
                $end  =$cols_1[3]+1;
            }else{
                print STDERR "@cols_1\n";
                print STDERR "@cols_2\n";
                next;
                #exit;
            }

            # remove PCR duplicates
            if( $remove_duplicate == 1 ){
                if( not defined $href_check_duplicate->{$cols_1[2]}->{"$start,$end"} ){
                    $href_check_duplicate->{$cols_1[2]}->{"$start,$end"}=1;
                }else{
                    $href_check_duplicate->{$cols_1[2]}->{"$start,$end"}++;
                    next;
                }
            }

            # create index for start and end, realize fast query
            my $idx_start=int($start/1000)*1000;
            my $idx_end  =(int($end/1000)+1)*1000;

            push @{$href_frag->{$cols_1[2]}->{"$idx_start\t$idx_end"}},[$start,$end];

        }else{
                #print STDERR "@cols_1\n";
                #print STDERR "@cols_2\n";
                next;
                #exit;
        }

    } 

    close IN;
    return $href_frag;
}
