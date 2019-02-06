#!/bin/sh
BINDIR=`dirname "${BASH_SOURCE[0]}"`
BINDIR=`readlink -f ${BINDIR}`
echo "BINDIR: " $BINDIR

MAINDIR=`pwd`
MAINDIR="$MAINDIR/process_DSB-seq_data"
echo "MAINDIR: " $MAINDIR
if [ ! -d $MAINDIR ]; then
    mkdir $MAINDIR
fi

EXPECTED_ARGS=6

if [ $# -lt ${EXPECTED_ARGS} ]
then
  echo "
Usage: `basename $0` <reference name> <chromosome length> <reference bowtie index> <enzyme name> <enzyme type: 5, 3, or blunt> <enzyme cutting sites> <DSB sequencing reads> <output prefix>
  "
  exit
fi

reference_name=$1 
chr_length=`readlink -f $2`
bowtie_index=`readlink -f $3 `
enzyme_name=$4
enzyme_type=$5
enzyme_cutting_sites=`readlink -f $6`
reads=`readlink -f $7`
Prefix=$8

##################### set which step you want to run #####################

# it will map DSB sequencing reads to genome
map_genome=1
# it will create wig file from btt file
create_wig=1
# count reads from enzyme cutting sites
count_enzyme_reads=1

##################### run scripts #####################

# map to genome
btt_dir=map_SE_to_$reference_name\_btt
if [ $map_genome == 1 ]
then
        echo ""
        echo "Running: map to genome, create btt file, please don't terminate the running software"
        echo ""
        cd $MAINDIR
        #echo sh $BINDIR/others/mapping_SE.sh $reference_name $bowtie_index \"-m1 -v1 -5 0 -3 0 -p2 -r\" $btt_dir $Prefix close_barcode $reads 
        #echo ""
        sh $BINDIR/others/mapping_SE.sh $reference_name $bowtie_index "-m1 -v1 -5 0 -3 0 -p2 -r" $btt_dir $Prefix close_barcode $reads > $btt_dir.log 2>&1
 
        echo "Done"

fi

# create wig 
depth_dir=$btt_dir\_depth
if [ $create_wig == 1 ]
then
        echo ""
        echo "Running: convert genome btt to wig"
        echo ""
        cd $MAINDIR
        if [ ! -d $depth_dir ]; then
            mkdir $depth_dir
        fi
        cd $depth_dir

        # create startDepth.wig, startDepth.plus_and_minus.wig, startDepth.plus-minus.wig
        #echo $BINDIR/PERL/btt_to_depth_fast.pl ../$btt_dir/$Prefix\.bt.btt $chr_length $Prefix 1 1 
        #echo ""
        perl $BINDIR/PERL/btt_to_depth_fast.pl ../$btt_dir/$Prefix\.bt.btt $chr_length $Prefix 1 1 > $Prefix\.btt_to_depth.log 2>&1

        echo "Done"
fi

# count enzyme reads 
enzyme_reads_dir=get_$enzyme_name\_reads
if [ $count_enzyme_reads == 1 ]
then
        echo ""
        echo "Running: count reads from $enzyme_name enzyme cutting sites"
        echo ""
        cd $MAINDIR

        #echo $BINDIR/PERL/countReads.pl $enzyme_cutting_sites $depth_dir/$Prefix.startDepth.plus_and_minus.wig -d $enzyme_reads_dir -p $Prefix.$enzyme_name 
        #echo ""
        perl $BINDIR/PERL/countReads.pl $enzyme_cutting_sites $depth_dir/$Prefix.startDepth.plus_and_minus.wig -d $enzyme_reads_dir -p $Prefix.$enzyme_name > $enzyme_reads_dir.log 2>&1

        echo "Done"
fi

