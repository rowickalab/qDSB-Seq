#!/usr/bin/sh

BINDIR=`dirname "${BASH_SOURCE[0]}"`

EXPECTED_ARGS=6

if [ $# -lt ${EXPECTED_ARGS} ]
then
  echo "
Usage: `basename $0` <organism> <bowtie_index> <bowtie parameters> <outdir> <bt_prefix> <hits_directory> <Read1.sort.fasta> <Read2.sort.fasta>

<organism>        human, mouse, yeast, yeast_qDSB-seq 
<bowtie_index>    bowtie index 
                  for example:
                  yeast qDSB-seq: /data/store/yizhu/my_databases/yeast/sc3_generate_G_MEC/sc3_generated_G_MEC.bowtie
                  yeast sc3: /data/store/yizhu/my_databases/yeast/sc2011/sc3.bowtie
                  human: /data/store/yizhu/my_databases/human/bowtie_index/hg19
                  mouse: /data/store/yizhu/my_databases/mouse/mm10/bowtie_index/mm10 
<bowtie_parameters>             bowtie parameters: '-m1 -v1 -5 0 -3 10 -p5 -[fqr]'
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

BOWTIE=bowtie
BTT=$BINDIR/../src/btt

OUTPUTBT=$OUTDIR/${PREFIX}.bt
OUTPUTBTT=$OUTDIR/${PREFIX}.bt.btt
OUTDIRHITS=$OUTDIR/${HITSDIR}

mkdir ${OUTDIR}
echo "$OUTDIR has been created!"
mkdir ${OUTDIRHITS}
echo "$OUTDIRHITS has been created!"

echo "Running Bowtie"

if [[ -f $FILE1 && -f $FILE2 ]]
then

    $BOWTIE $BOWTIE_INDEX $BOWTIE_PARAMETERS $FILE1,$FILE2 > $OUTPUTBT 2>$OUTPUTBT.log 

elif [[ -f $FILE1 ]] 
then

    $BOWTIE $BOWTIE_INDEX $BOWTIE_PARAMETERS --un $OUTPUTBT\.unmap.fa $TYPE $FILE1 > $OUTPUTBT 2>$OUTPUTBT.log #or -v 1

else

    echo "${FILE1} is not exist!"
    exit

fi
    
if [[ -f ${OUTPUTBT} ]]
then

    echo "Convert bt to btt using btt" 
    ${BTT} ${OUTPUTBT}

fi

#check btt
if [ ! -f ${OUTPUTBTT} ]
then
    echo "${OUTPUTBTT} is not exist!"
    exit
fi

