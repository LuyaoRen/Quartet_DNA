#!/bin/bash

mkdir tardis
cd tardis
sample=LCL5
mkdir ${sample}
cd ${sample}
ln -s ${sample}_bwa_chr1_22_XYM.bam
ln -s ${sample}_bwa_chr1_22_XYM.bam.bai

map=bwa
callSV=tardis

fasta=chr1_22_XYM.fa
GRCh38_sonic=$path/GRCh38.sonic

tardis-nocram -i ${sample}_bwa_chr1_22_XYM.bam --ref $fasta --sonic ${GRCh38_sonic}  \
--out ${sample}_${map}_${callSV}