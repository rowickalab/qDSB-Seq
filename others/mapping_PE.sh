#!/usr/bin/sh

BINDIR=`dirname "${BASH_SOURCE[0]}"`
BINDIR=`readlink -f ${BINDIR}`
BOWTIE=bowtie
BTT=$BINDIR/../src/btt

EXPECTED_ARGS=6

if [ $# -lt ${EXPECTED_ARGS} ]
then
  echo "
Usage: `basename $0` <organism> <bowtie_index> <bowtie parameters> <outdir> <bt_prefix> <hits_directory> <Read1.sort.fasta> <Read2.sort.fasta>

<organism>        human, mouse, yeast
<bowtie_index>    bowtie index
                  for example:
                  yeast qDSB-seq paper: /data/store/yizhu/my_databases/yeast/sc3_generate_G_MEC/sc3_generated_G_MEC.bowtie
                  yeast sc3: /data/store/yizhu/my_databases/yeast/sc2011/sc3.bowtie
                  human: /data/store/yizhu/my_databases/human/bowtie_index/hg19
                  human38: /data/store/yizhu/my_databases/human/hg38/hg38.bowtie
                  mouse: /data/store/yizhu/my_databases/mouse/mm10/bowtie_index/mm10 
                  yeast rDNA: /data/store/yizhu/my_databases/yeast/rDNA/from_Benjamin/rDNA_9137bp_dig_BigIII.bowtie
<bowtie_parameters>             bowtie parameters: '-m1 -v1 -5 0 -3 20 -p5 --maxins 1000 -[fqr]'
<outdir>          The output directory of all outputs
<bt_prefix>       The output prefix of .bt,.btt
<hits_directory>  The directory of hits files: no_barcode, close_barcode,distant_barcode
<Read1>           Read 1 file, if it's single-end mapping, you can use multiple files and concatenate by ',' 
<Read2>           Read 2 file for paired-end mapping

"
  exit
fi

ORGANISM=$1
BOWTIE_INDEX=$2
BOWTIE_PARAMETERS=$3
OUTDIR=$4
PREFIX=$5
HITSDIR=$6
FILE1=$7
FILE2=$8

echo $ORGANISM
echo $BOWTIE_INDEX
echo $BOWTIE_PARAMETERS
echo $OUTDIR
echo $PREFIX
echo $HITSDIR
echo $FILE1
echo $FILE2

OUTPUTBT=$OUTDIR/${PREFIX}.bt
OUTPUTBTT=$OUTDIR/${PREFIX}.bt.btt
OUTDIRHITS=$OUTDIR/${HITSDIR}

mkdir ${OUTDIR}
mkdir ${OUTDIRHITS}

echo "Running Bowtie"

if [[ -f $FILE1 && -f $FILE2 ]]
then
    echo "$BOWTIE $BOWTIE_INDEX $BOWTIE_PARAMETERS -1 $FILE1 -2 $FILE2 > $OUTPUTBT 2>$OUTPUTBT.log "
    $BOWTIE $BOWTIE_INDEX $BOWTIE_PARAMETERS -1 $FILE1 -2 $FILE2 > $OUTPUTBT 2>$OUTPUTBT.log 
else
    echo "${FILE1} or ${FILE2} is not exist!"
    exit
fi
    
if [[ -f ${OUTPUTBT} ]]
then
    echo "Convert bt to btt using btt" 
    ${BTT} ${OUTPUTBT}
fi
exit

#check btt
if [ ! -f ${OUTPUTBTT} ]
then
    echo "${OUTPUTBTT} is not exist!"
    exit
fi

exit

