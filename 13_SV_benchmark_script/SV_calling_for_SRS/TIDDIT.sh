#!/bin/bash

mkdir TIDDIT
cd TIDDIT

bam_1=LCL5_bwa_chr1_22_XYM.bam
bam_2=LCL6_bwa_chr1_22_XYM.bam
bam_3=LCL7_bwa_chr1_22_XYM.bam
bam_4=LCL8_bwa_chr1_22_XYM.bam

sample=LCL5
mkdir ${sample}
cd ${sample}
threads=2
map=bwa
callSV=TIDDIT
fasta=chr1_22_XYM.fa

source activate py37
TIDDIT.py --sv --bam ${bam_1} --ref ${fasta} -o ${sample}_${map}_${callSV} -n 2 --force_ploidy 
