#!/bin/bash

cd 
mkdir manta
cd manta

bam_1=Quartet_DNA_ILM_Nova_ARD_LCL5_1_20181108.sorted.deduped.bam
bam_2=Quartet_DNA_ILM_Nova_ARD_LCL6_1_20181108.sorted.deduped.bam
bam_3=Quartet_DNA_ILM_Nova_ARD_LCL7_1_20181108.sorted.deduped.bam
bam_4=Quartet_DNA_ILM_Nova_ARD_LCL8_1_20181108.sorted.deduped.bam

sample=LCL5
mkdir ${sample}
cd ${sample}
threads=2

samtools view -h ${bam_1}|grep -E -v "chrUn|random|SV40|MCV|HTLV|KSHV|HIV|HCV|HBV|CMV|HPV|EBV" -|samtools view -bhS -> ${sample}_bwa_chr1_22_XYM.bam
samtools index -@ $threads ${sample}_bwa_chr1_22_XYM.bam

fasta=chr1_22_XYM.fa
outdir=manta/${sample}

configManta.py --bam=${sample}_bwa_chr1_22_XYM.bam --referenceFasta=${fasta} --runDir=${outdir}
manta/${sample}/runWorkflow.py -j 4 -g 10

