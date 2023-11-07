#SV merging for each sample
ls *LCL5*_sorted.vcf >> LCL5_merge_pacbio12_ont.txt
ls *LCL6*_sorted.vcf >> LCL6_merge_pacbio12_ont.txt
ls *LCL7*_sorted.vcf >> LCL7_merge_pacbio12_ont.txt
ls *LCL8*_sorted.vcf >> LCL8_merge_pacbio12_ont.txt
for  i in `ls LCL*_merge_pacbio12_ont.txt`;do \
/Jasmine-1.0.1/run.sh \
max_dist=1000 min_support=0 threads=2 --allow_intrasample --keep_var_ids --ignore_strand  file_list=${i} out_file=${i%.*}_Intrasample.vcf;done

#2TechOr6Pipe
vcf=${sample}_merge_pacbio12_ont_Intrasample.vcf
perl -F'\t' -lane 'if(/#/){print}else{($id)=$F[7]=~m/IDLIST=([^\s]+)/;my %count;my @tech=();my @pipe=();@multiID=split/,/,$id;foreach $sample (@multiID){my @ID=split /_LCL/,$sample;push(@tech,$ID[0]);$MID=$sample;$MID=~s/_LCL[0-9]+//g;$MID=~s/_[0-9]+_[A-Z]+//g;push(@pipe,$MID)}if(scalar(grep{++$count{$_}<2} @tech)>=2 || scalar(grep{++$count{$_}<2} @pipe)>=6){print}else{next}}' \
$vcf>${vcf%.*}_2TechOr6Pipe.vcf

#SV merging for quartet
ls *_2TechOr6Pipe.vcf >> quartet.txt
/Jasmine-1.0.1/run.sh \
max_dist=1000 min_support=0 threads=2 --keep_var_ids --ignore_strand --output_genotypes file_list=quartet.txt out_file=quartet.vcf

#Removing SVs in centromere or gap
centro=centro_pericentro_region.bed
gap=hg38_gaps.bed
vcf=quartet.vcf
cat <(grep "#" ${vcf}) <(bedtools intersect -a ${vcf} -b <(cat $centro $gap|sort -k1,1 -k2,2n|cut -f 1-3|bedtools merge -i -) -v) \
> ${vcf%.*}_filter_centromere_gap.vcf

#less10M
perl -lane '{($LEN)=$F[7]=~m/SVLEN=([^;]+)/;if(abs($1)>=10000000){next}else{print}}' quartet_filter_centromere_gap.vcf \
> quartet_filter_centromere_gap_less10M.vcf

#INS/DEL
grep -E "#|SVTYPE=INS|SVTYPE=DEL" quartet_filter_centromere_gap_less10M.vcf\
>quartet_filter_centromere_gap_less10M_INS_DEL.vcf

#Getting read names of supporting SVs and allele sequences
#PB2_minimap2_bam
vcf=quartet_filter_centromere_gap_less10M_INS_DEL.vcf
grep -E SUPP_VEC=1 $vcf|sort -k1,1 -k2,2n ->quartet_ReadName_LCL5.vcf
grep -E SUPP_VEC=01 $vcf|sort -k1,1 -k2,2n ->quartet_ReadName_LCL6.vcf
grep -E SUPP_VEC=001 $vcf|sort -k1,1 -k2,2n ->quartet_ReadName_LCL7.vcf
grep -E SUPP_VEC=0001 $vcf|sort -k1,1 -k2,2n ->quartet_ReadName_LCL8.vcf
sniffles -t $threads -m ${PB2_bam} -n -1 -s 2 -v ${vcf%.*}_sniffles.vcf --Ivcf $vcf

#Refining allele sequences
java -cp /Iris/src Iris \
threads=2 \
max_ins_dist=100 \
max_len_change=0.25 \
max_out_length=1000000 \
--also_deletions\
--pacbio \
genome_in=$fasta vcf_in=$vcf reads_in=${PB2_bam} vcf_out=${vcf%.*}_Iris.vcf

#SV genotyping
cat *_Iris.vcf| \
perl -lane 'if(/SVTYPE=DEL/){if(length($F[3])>=50){print}else{next}}elsif(/SVTYPE=INS/){if(length($F[4])>=50){print}else{next}}else{next}' -|\
sort -k1,1 -k2,2n|awk -v OFS="\t" '{print $1,$2,$3,toupper($4),toupper($5),".","PASS",$8}' -|\
grep -v "NNNNNN" -|grep -E -v "<INS>|<DEL>" -|perl -lane 's/;RNAMES=[^;]+//g;print' -\
>quartet_filter_centromere_gap_less10M_INS_DEL_sniffles_Iris_Len50.vcf
LRcaller -fa ${fasta} -ora -nt $threads -a seqan ${bam} ${sv} ${sv_geno}
sniffles -t $threads -m ${bam} -v ${sv_geno} --Ivcf ${sv}
python3 svjedi.py -t $Threads -d pb -ms 3 -v $sv -r $ref -i $fasta -o ${sv_geno}
python3 svjedi.py -t $Threads -d ont -ms 3 -v $sv -r $ref -i $fasta -o ${sv_geno}
/Jasmine-1.0.1/run.sh \
max_dist=1000 min_support=0 threads=2 --keep_var_ids --ignore_strand --output_genotypes \
file_list=${sample}_svjedi_sniffles_LRcaller.txt out_file=${sample}_svjedi_sniffles_LRcaller.vcf
perl -lane 'my $num=6;if(/#/){print}else{if(/(0\/0:.*){$num,}/){s/\S\/\S/0\/0/g;print $_}elsif(/(1\/1:.*){$num,}/){s/\S\/\S/1\/1/g;print $_}elsif(/(0\/1:.*){$num,}/){s/\S\/\S/0\/1/g;print $_}}' \
$vcf|cut -f 1-10 ->${vcf%.*}_Geno6.vcf
/Jasmine-1.0.1/run.sh \
max_dist=1000 min_support=0 threads=2 --keep_var_ids --ignore_strand --output_genotypes \
file_list=quartet_geno.txt out_file=quartet_geno.vcf

#SV benchmark set
header=Jasmine_header_mfc.txt
vcf=quartet_geno.vcf
bcftools +mendelian <(cat $header <(grep -E "SUPP_VEC=1111" $vcf|grep -E -v "\.\/\."|grep -E -v "0\/0:.*0\/0:.*0\/0:.*0\/0:.*"| \
grep -v -E "chrX|chrY|chrM"| \
awk -v OFS="\t" '/#/{print}{split($10,LCL5,":");split($11,LCL6,":")}{if(LCL5[1]==LCL6[1]){print $0}}' -| \
perl -lane 'if(/^#/){print}else{print join("\t",@F[0..8],$F[12],$F[11],$F[9])}')) \
 -t mother,father,child -l+ > quartet_geno_Mendelian.vcf