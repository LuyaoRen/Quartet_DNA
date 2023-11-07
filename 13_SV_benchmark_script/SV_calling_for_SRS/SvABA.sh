#!/bin/bash

mkdir SvABA
cd SvABA

bam_1=LCL5_bwa_chr1_22_XYM.bam
bam_2=LCL6_bwa_chr1_22_XYM.bam
bam_3=LCL7_bwa_chr1_22_XYM.bam
bam_4=LCL8_bwa_chr1_22_XYM.bam

sample=LCL5
mkdir ${sample}
cd ${sample}
threads=4
map=bwa
callSV=SvABA
fasta=chr1_22_XYM.fa

svaba run -t ${bam_1} -G ${fasta} -a ${sample} -p $threads --germline 
