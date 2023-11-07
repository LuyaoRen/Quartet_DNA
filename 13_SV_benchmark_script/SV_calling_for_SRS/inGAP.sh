#!/bin/bash

mkdir inGAP
cd inGAP

bam_1=LCL5_bwa_chr1_22_XYM.bam
bam_2=LCL6_bwa_chr1_22_XYM.bam
bam_3=LCL7_bwa_chr1_22_XYM.bam
bam_4=LCL8_bwa_chr1_22_XYM.bam

sample=LCL5
mkdir ${sample}
cd ${sample}
threads=2
map=bwa
callSV=inGAP
fasta=chr1_22_XYM.fa

chrom=chr$chrnum

samtools view ${sample}_bwa_SVseq2_${chrom}.bam>${sample}_bwa_SVseq2_${chrom}.sam
java -mx20000m -jar /inGAP_3_1_1/inGAP.jar SVP -SE 5 -PE 5 -SIZE 1000000 -r ${fasta} -i ${sample}_bwa_SVseq2_${chrom}.sam -o ${sample}_${map}_${callSV}_${chrom}.txt
