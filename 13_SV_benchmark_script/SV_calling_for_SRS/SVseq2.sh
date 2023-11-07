#!/bin/bash

mkdir SVseq2
cd SVseq2

bam_1=LCL5_bwa_chr1_22_XYM.bam
bam_2=LCL6_bwa_chr1_22_XYM.bam
bam_3=LCL7_bwa_chr1_22_XYM.bam
bam_4=LCL8_bwa_chr1_22_XYM.bam

sample=LCL5
mkdir ${sample}
cd ${sample}
threads=2
map=bwa
callSV=SVseq2

chrom=chr$chrnum
samtools view -@ ${threads} -bh -F 256 -f 2 ${bam_1} ${chrom} > ${sample}_${map}_${callSV}_${chrom}.bam
samtools index ${sample}_${map}_${callSV}_${chrom}.bam
chr=chroms/${chrom}.fa

#Call DELs
SVseq2_2 -r $chr -c ${chrom} -b ${sample}_${map}_${callSV}_${chrom}.bam --o ${sample}_${map}_${callSV}_${chrom}_DEL.txt --is insert std

#Call INSs 
SVseq2_2 -insertion -c ${chrom} -b ${sample}_${map}_${callSV}_${chrom}.bam --o ${sample}_${map}_${callSV}_${chrom}_INS.txt
