#!/bin/bash

mkdir wham
cd wham

bam_1=LCL5_bwa_chr1_22_XYM.bam
bam_2=LCL6_bwa_chr1_22_XYM.bam
bam_3=LCL7_bwa_chr1_22_XYM.bam
bam_4=LCL8_bwa_chr1_22_XYM.bam

sample=LCL5
mkdir ${sample}
cd ${sample}
threads=4
map=bwa
callSV=wham
fasta=chr1_22_XYM.fa

WHAM-GRAPHENING -f ${bam_2} -a ${fasta} -x ${threads} -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM > ${sample}_${map}_${callSV}.vcf
