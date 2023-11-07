
#Read mapping
name=quartet_subreads_name

#minimap2
minimap2 -ax map-pb --MD -Y -t $threads $reference ${name}.fa | samtools view -bhS -@ $threads -o ${name}.bam
samtools sort -@ $threads -O BAM -o ${name}_sorted.bam ${name}.bam
samtools index -@ $threads ${name}_sorted.bam

#ngmlr
ngmlr -t $threads -r $reference -q ${name}.fa -x pacbio | samtools view -bhS -@ $threads - > ${name}.bam
samtools sort -@ $threads -O BAM -o ${name}_sorted.bam ${name}.bam
samtools index -@ $threads ${name}_sorted.bam

#pbmm2
pbmm2 align $reference ${sample_dir}/$name ${name}.bam --sort --median-filter --sample $sample
samtools index -@ $threads ${name}.bam
