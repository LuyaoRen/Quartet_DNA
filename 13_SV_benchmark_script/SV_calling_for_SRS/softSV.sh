#!/bin/bash

mkdir SoftSV
cd SoftSV

bam_1=LCL5_bwa_chr1_22_XYM.bam
bam_2=LCL6_bwa_chr1_22_XYM.bam
bam_3=LCL7_bwa_chr1_22_XYM.bam
bam_4=LCL8_bwa_chr1_22_XYM.bam

sample=LCL5
mkdir ${sample}
cd ${sample}
threads=2
map=bwa
callSV=SoftSV
fasta=chr1_22_XYM.fa


chrom=chr$chrnum
mkdir ${sample}_${map}_${callSV}
samtools view -@ ${threads} -bh ${bam_1} ${chrom} > ${sample}_${map}_${callSV}_${chrom}.bam
samtools index ${sample}_${map}_${callSV}_${chrom}.bam
SoftSV -i ${sample}_${map}_${callSV}_${chrom}.bam -r ${chrom} -o ${sample}_${map}_${callSV} -p ${sample}_${map}_${callSV}_${chrom}
