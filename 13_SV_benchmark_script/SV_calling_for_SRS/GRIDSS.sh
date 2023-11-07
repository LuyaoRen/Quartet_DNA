#!/bin/bash

mkdir GRIDSS
cd GRIDSS

bam_1=LCL5_bwa_chr1_22_XYM.bam
bam_2=LCL6_bwa_chr1_22_XYM.bam
bam_3=LCL7_bwa_chr1_22_XYM.bam
bam_4=LCL8_bwa_chr1_22_XYM.bam

sample=LCL5
mkdir ${sample}
cd ${sample}
threads=4
map=bwa
callSV=GRIDSS

fasta=chr1_22_XYM.fa
outdir=GRIDSS/${sample}

gridss.sh \
--reference $fasta --output ${sample}_${map}_${callSV}.vcf.gz \
--assembly ${sample}_${map}_${callSV}.bam \
--threads $threads \
--jar /software/gridss-2.8.1-gridss-jar-with-dependencies.jar \
--workingdir $outdir \
--jvmheap 25g \
--blacklist GRIDSS/hg38_gaps_chr1_22_XYM.bed \
--steps All \
--maxcoverage 50000 \
--labels ${sample} ${bam_1}
