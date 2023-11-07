#!/bin/bash

threads=2
source activate py27
cd 
mkdir lumpy
cd lumpy
sample=LCL5
mkdir ${sample}
cd ${sample}
ln -s ${sample}_bwa_chr1_22_XYM.bam

#samtools_v0.1.19
# Extract the discordant paired-end alignments.
samtools_v0.1.19 view -@ $threads -b -F 1294 ${sample}_bwa_chr1_22_XYM.bam > ${sample}.discordants.unsorted.bam

# Extract the split-read alignments
samtools_v0.1.19 view -@ $threads -h ${sample}_bwa_chr1_22_XYM.bam \
    | /public/home/fan_lab/xkduan/software/lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin \
    | samtools_v0.1.19 view -@ $threads -Sb - \
    > ${sample}.splitters.unsorted.bam

# Sort both alignments
samtools_v0.1.19 sort -@ $threads ${sample}.discordants.unsorted.bam ${sample}.discordants
samtools_v0.1.19 sort -@ $threads ${sample}.splitters.unsorted.bam ${sample}.splitters


/public/home/fan_lab/xkduan/software/lumpy-sv/bin/lumpyexpress \
    -B ${sample}_bwa_chr1_22_XYM.bam \
    -S ${sample}.splitters.bam \
    -D ${sample}.discordants.bam \
    -o ${sample}_lumpy.vcf
