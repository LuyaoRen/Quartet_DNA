
#SV calling
#cuteSV
#minimap2/ngmlr
cuteSV --threads ${threads} --genotype -s 3 ${bam} ${sv} ${work_dir}
#cutesv_3var
for i in `ls *_cutesv.vcf`;do grep -v -e "DEL/INV" -e "DUP/INS" -e "INVDUP" $i| \
perl -lane 'if(/SVTYPE=INS;.*?SVLEN=([-|\d]+)/ || /SVTYPE=DEL;.*?SVLEN=([-|\d]+)/ || /SVTYPE=DUP;.*?SVLEN=([-|\d]+)/ || /SVTYPE=INV;.*?SVLEN=([-|\d]+)/){if(abs($1)>=50){print}}else{print}' -| \
reset_var_depth.sh sniffles - ${i%.*}_50bp_5sv_3var.vcf;grep -v "#" -c ${i%.*}_50bp_5sv_3var.vcf;done

#sniffles
#minimap2/ngmlr
sniffles -t $threads -s 3 -m ${bam} -v ${sv}
#pbmm2
samtools calmd --output-fmt bam -@ $threads ${PB2_bam} ${fasta} > ${PB2_bam%.*}_MD.bam
sniffles -t $threads -s 2 -m ${PB2_bam%.*}_MD.bam -v ${sv}
#sniffles_3var
for i in `ls *_sniffles.vcf`;do grep -v -e "DEL/INV" -e "DUP/INS" -e "INVDUP" $i| \
perl -lane 'if(/SVTYPE=INS;.*?SVLEN=([-|\d]+)/ || /SVTYPE=DEL;.*?SVLEN=([-|\d]+)/ || /SVTYPE=DUP;.*?SVLEN=([-|\d]+)/ || /SVTYPE=INV;.*?SVLEN=([-|\d]+)/){if(abs($1)>=50){print}}else{print}' -| \
reset_var_depth.sh sniffles - ${i%.*}_50bp_5sv_3var.vcf;grep -v "#" -c ${i%.*}_50bp_5sv_3var.vcf;done

#nanosv
#minimap2/ngmlr
NanoSV --bed hg38_genome_sample_${chrom}.bed -t ${threads} -s samtools ${sample}_${map}_nanosv_${chrom}.bam -o ${sample}_${map}_nanosv_${chrom}.vcf
for i in `ls *_nanosv*.vcf`;do perl -lane 'if(/SVTYPE=(.*?);/){print $1}' $i|sort|uniq -c;\
perl -lane 's/;ALT_READ_IDS=[^\s]+//g;print' $i>>../${sample}_${map}_nanosv.vcf;done
#nanosv_3var
for i in `ls *vcf`;do grep -c -v "#" $i; \
grep -v "#" $i|perl -lane 's/;ALT_READ_IDS=[^\s]+//g;print' -|\
perl -lane 'if(/SVTYPE=INS;.*?SVLEN=([-|\d]+)/ || /SVTYPE=DEL;.*?SVLEN=([-|\d]+)/ || /SVTYPE=DUP;.*?SVLEN=([-|\d]+)/ || /SVTYPE=INV;.*?SVLEN=([-|\d]+)/){if(abs($1)>=50){print}}else{print}' -|\
perl -lane 'my $var=3;if(/^#/){print $_}else{@format=split /:/,$F[9];@allele=split /,/,$format[2];if($allele[1]>=$var){print $_}}' -\
> ${i%.*}_50bp_5sv_3var.vcf;done

#svim
#minimap2/ngmlr
for i in `ls *_svim.vcf`;do sed -e 's/DUP_INT/DUP/g' -e 's/DUP:TANDEM/DUP/g' $i| \
perl -lane 'if(/SVTYPE=INS;.*?SVLEN=([-|\d]+)/ || /SVTYPE=DEL;.*?SVLEN=([-|\d]+)/ || /SVTYPE=DUP;.*?SVLEN=([-|\d]+)/){if(abs($1)>=50){print}}elsif(/SVTYPE=INV;.*?END=([-|\d]+)/){if(abs($1-$F[1])>=50){print}}else{print}' -| \
reset_var_depth.sh svim - ${i%.*}_50bp_5sv_3var.vcf;grep -v "#" -c ${i%.*}_50bp_5sv_3var.vcf;done

#pbsv
#pbmm2/ngmlr
#TAG='@RG\tID:LCL5_ngmlr_m\tSM:LCL5\tLB:ngmlr'
#samtools addreplacerg -@ $threads -m overwrite_all -r $TAG ./${sample}_ngmlr.bam -O BAM > ${sample}_${map}_rg.bam
pbsv discover --tandem-repeats $trf ${sample}_${map}.bam ${sample}_${map}.svsig.gz
pbsv call ${fasta} ${sample}_${map}.svsig.gz ${sample}_${map}_pbsv.vcf

for i in `ls *_pbsv.vcf`;do grep -v "cnv" $i| \
perl -lane 'if(/SVTYPE=INS;.*?SVLEN=([-|\d]+)/ || /SVTYPE=DEL;.*?SVLEN=([-|\d]+)/ || /SVTYPE=DUP;.*?SVLEN=([-|\d]+)/){if(abs($1)>=50){print}}elsif(/SVTYPE=INV;.*?END=([-|\d]+)/){if(abs($1-$F[1])>=50){print}}else{print}' -| \
reset_var_depth.sh pbsv - ${i%.*}_50bp_5sv_3var.vcf;grep -v "#" -c ${i%.*}_50bp_5sv_3var.vcf;done