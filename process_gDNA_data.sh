#!/bin/sh
BINDIR=`dirname "${BASH_SOURCE[0]}"`
BINDIR=`readlink -f ${BINDIR}`
echo "BINDIR: " $BINDIR

MAINDIR=`pwd`
MAINDIR="$MAINDIR/process_gDNA_data"
echo "MAINDIR: " $MAINDIR
if [ ! -d $MAINDIR ]; then
    mkdir $MAINDIR
fi

EXPECTED_ARGS=10

if [ $# -lt ${EXPECTED_ARGS} ]
then
  echo "
Usage: `basename $0` <reference name> <chromosome length> <reference bowtie index> <enzyme name> <enzyme type: 5, 3, or blunt> <enzyme cutting sites> <gDNA reads R1> <gDNA reads R2> <background_list> <output prefix>
  "
  exit
fi

reference_name=$1 
chr_length=`readlink -f $2`
bowtie_index=`readlink -f $3 `
enzyme_name=$4
enzyme_type=$5
enzyme_cutting_sites=`readlink -f $6`
R1=`readlink -f $7`
R2=`readlink -f $8`
background_list=`readlink -f $9`
Prefix=${10}

##################### set which step you want to run #####################

# it will map DSB sequencing reads to genome
map_genome=1
# compute cutting efficiency
compute_cutting_efficiency=1

##################### run scripts #####################

# map to genome
btt_dir=map_PE_to_$reference_name\_btt
if [ $map_genome == 1 ]
then
        echo ""
        echo "Running: map to genome, create btt file, please don't terminate the running software"
        echo ""
        cd $MAINDIR

        #echo sh $BINDIR/others/mapping_PE.sh $reference_name $bowtie_index \"-m1 -v1 -5 0 -3 0 -p2 --maxins 1000 -r\" $btt_dir $Prefix no_barcode $R1 $R2
        #echo ""
        sh $BINDIR/others/mapping_PE.sh $reference_name $bowtie_index "-m1 -v1 -5 0 -3 0 -p2 --maxins 1000 -r" $btt_dir $Prefix no_barcode $R1 $R2 >$btt_dir.log 2>&1

        perl $BINDIR/PERL/split_btt.pl $btt_dir/$Prefix.bt.btt 

        echo "Done"

fi

# compute cutting efficiency
fcut_dir=cutting_efficiency_$enzyme_name
if [ $compute_cutting_efficiency == 1 ]
then
        echo ""
        echo "Running: compute $enzyme_name cutting efficiency"
        echo ""
        cd $MAINDIR
        if [ ! -d $fcut_dir ]; then
            mkdir $fcut_dir
        fi
        cd $fcut_dir


        if [ $enzyme_type == 5 ]
        then
            #echo  perl $BINDIR/PERL/compute_cutting_efficiency.pl $enzyme_cutting_sites ../$btt_dir/$Prefix.bt.btt -l 0 -r 0 -exact -switch -prefix $Prefix.$enzyme_name
            #echo ""
            perl $BINDIR/PERL/compute_cutting_efficiency.pl $enzyme_cutting_sites ../$btt_dir/$Prefix.bt.btt -l 0 -r 0 -exact -switch -prefix $Prefix.$enzyme_name > $Prefix.$enzyme_name.log 2>&1

            perl $BINDIR/PERL/compute_cutting_efficiency.pl $background_list ../$btt_dir/$Prefix.bt.btt -l 0 -r 3 -exact -switch -prefix $Prefix.$enzyme_name.background > $Prefix.$enzyme_name.background.log 2>&1
            mv $Prefix.$enzyme_name.background.*.exact.switch.txt $Prefix.$enzyme_name.background.exact.switch.txt
        else
            #echo  perl $BINDIR/PERL/compute_cutting_efficiency.pl $enzyme_cutting_sites ../$btt_dir/$Prefix.bt.btt -l 0 -r 0 -exact -prefix $Prefix.$enzyme_name
            #echo ""
            perl $BINDIR/PERL/compute_cutting_efficiency.pl $enzyme_cutting_sites ../$btt_dir/$Prefix.bt.btt -l 0 -r 0 -exact -prefix $Prefix.$enzyme_name > $Prefix.$enzyme_name.log 2>&1

            perl $BINDIR/PERL/compute_cutting_efficiency.pl $background_list ../$btt_dir/$Prefix.bt.btt -l 0 -r 3 -exact -prefix $Prefix.$enzyme_name.background > $Prefix.$enzyme_name.background.log 2>&1
            mv $Prefix.$enzyme_name.background.*.exact.switch.txt $Prefix.$enzyme_name.background.exact.switch.txt
            
        fi

        echo "Done"
fi
