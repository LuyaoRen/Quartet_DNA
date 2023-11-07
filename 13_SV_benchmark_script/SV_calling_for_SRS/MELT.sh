#!/bin/bash

mkdir MELT
cd MELT

bam_1=LCL5_bwa_chr1_22_XYM.bam
bam_2=LCL6_bwa_chr1_22_XYM.bam
bam_3=LCL7_bwa_chr1_22_XYM.bam
bam_4=LCL8_bwa_chr1_22_XYM.bam

sample=LCL5
mkdir ${sample}
cd ${sample}
threads=2
map=bwa
callSV=MELT
fasta=chr1_22_XYM.fa

gene=/MELTv2.2.0/add_bed_files/Hg38/Hg38.genes.bed
ALU=/MELTv2.2.0/me_refs/Hg38/ALU_MELT.zip
SVA=/MELTv2.2.0/me_refs/Hg38/SVA_MELT.zip
LINE1=/MELTv2.2.0/me_refs/Hg38/LINE1_MELT.zip
HERVK=/MELTv2.2.0/me_refs/Hg38/HERVK_MELT.zip


#Preprocess
DIR=/MELTv2.2.0
java -Xmx2G -jar $DIR/MELT.jar Preprocess -bamfile ${bam_1} -h ${fasta}

#Call SVs 
java -Xmx2G -jar $DIR/MELT.jar Single -bamfile ${bam_1} -h ${fasta} -n $gene -t $ALU -w . -r 150 -e 420 -d 40000000 -c 30
java -Xmx2G -jar $DIR/MELT.jar Single -bamfile ${bam_1} -h ${fasta} -n $gene -t $SVA -w . -r 150 -e 420 -d 40000000 -c 30
java -Xmx2G -jar $DIR/MELT.jar Single -bamfile ${bam_1} -h ${fasta} -n $gene -t $LINE1 -w . -r 150 -e 420 -d 40000000 -c 30
java -Xmx2G -jar $DIR/MELT.jar Single -bamfile ${bam_1} -h ${fasta} -n $gene -t $HERVK -w . -r 150 -e 420 -d 40000000 -c 30
