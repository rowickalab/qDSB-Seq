perl ../qDSB-seq.pl test_i-BLESS.seq test_gDNA.R1.seq test_gDNA.R2.seq -s test -r G1 -g yeast -f reference_genome/test.reference_genome.fas -i reference_genome/test.reference_genome.bowtie -e NotI -t 5 -c NotI.bed -b background.bed -p output

#rm -r output.* process_*
