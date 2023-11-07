
#Confident regions
#Twins
grep "^R" D5_F7_asm.var.txt|awk '{sum+=$4-$3}END{print "Sum = ",sum}'
Sum =  2742426882
grep "^R" D5_M8_asm.var.txt|awk '{sum+=$4-$3}END{print "Sum = ",sum}'
Sum =  2740640768

cd /Quartet_data/SV_combine/pacbio_ont_pacbio2_combine_3/Assembly_SV_Region/paftools
#Find Diploid Assembled Regions: 
grep "^R" ./D5_haplotypeF7/D5_F7_asm.var.txt|cut -f 2-4 > D5_F7_one_contig.bed
grep "^R" ./D5_haplotypeM8/D5_M8_asm.var.txt|cut -f 2-4 > D5_M8_one_contig.bed
grep "^R" ./D6_haplotypeF7/D6_F7_asm.var.txt|cut -f 2-4 > D6_F7_one_contig.bed
grep "^R" ./D6_haplotypeM8/D6_M8_asm.var.txt|cut -f 2-4 > D6_M8_one_contig.bed
cat D5_F7_one_contig.bed D5_M8_one_contig.bed|sort -k1,1 -k2,2n|bedtools merge -i - > D5_F7M8_one_contig.bed
cat D6_F7_one_contig.bed D6_M8_one_contig.bed|sort -k1,1 -k2,2n|bedtools merge -i - > D6_F7M8_one_contig.bed


#Combine Regions:
bedtools intersect -a D5_F7M8_one_contig.bed -b D6_F7M8_one_contig.bed| \
grep -E -v "chrX|chrY" > D5_D6_F7M8_one_contig_autosomes.bed

#Discover Variants: & Compare Variants:
cd /Quartet_data/SV_combine/pacbio_ont_pacbio2_combine_3/Assembly_SV_Region/paftools/Confident_regions/
/software/Jasmine-1.0.1/run.sh  \
max_dist=1000 min_support=0 threads=2 --keep_var_ids --ignore_strand \
file_list=Tier1_paftools.txt out_file=Tier1_paftools.vcf
"#"|sort -k1,1 -k2,2n|cut -f 1-3|bedops --range -500:500 --everything - > paftools_NoMatchTier1_plus500bp.bed

#SVwiden
mkdir SVwiden_original
cd SVwiden_original
grep -E "SUPP_VEC=0" ../Tier1_paftools.vcf|grep -E -v "chrX|chrY"|sort -k1,1 -k2,2n > paftools_NoMatchTier1.vcf
/Pacbio_command/vcf2bed.sh paftools_NoMatchTier1.vcf| \
grep -v "#"|sort -k1,1 -k2,2n|cut -f 1-4 -> paftools_NoMatchTier1.bed

ref=/hg38_ref/chr1_22_XYM.fa/chr1_22_XYM.fa
cat <(grep "#" ../Tier1_paftools.vcf) paftools_NoMatchTier1.vcf > paftools_NoMatchTier1_h.vcf

svanalyzer widen --ref $ref --variants paftools_NoMatchTier1_h.vcf
perl -lane 'if(/REFWIDENED=(.*?):([0-9]+)-([0-9]+)/){print join("\t",$1,$2,$3,$F[2])}' widened.vcf| \
awk -v OFS='\t' '{if($3<=$2){print $1,$3,$2,$4}else{print}}' -| \
bedops --range -50:50 --everything - > widened_50bp.bed


#Trim Regions:
#NoPaftools
bed=/Quartet_data/SV_combine/pacbio_ont_pacbio2_combine_3/Assembly_SV_Region/paftools/D5_D6_F7M8_one_contig_autosomes.bed

bedtools subtract -a $bed -b <(cat paftools_NoMatchTier1.bed widened_50bp.bed|sort -k1,1 -k2,2n|bedops --range -50:50 --everything -)| \
sort -k1,1 -k2,2n|bedtools merge -i -|awk '{sum+=$3-$2}END{print "Sum =",sum}'
Sum = 2622728511

bedtools subtract -a $bed -b <(cat paftools_NoMatchTier1.bed widened_50bp.bed|sort -k1,1 -k2,2n|bedops --range -50:50 --everything -)| \
sort -k1,1 -k2,2n|bedtools merge -i -> Tier1_Confident_regions_autosomes_NoPaftools.bed


#F7/M8
cd /Quartet_data/SV_combine/pacbio_ont_pacbio2_combine_3/Assembly_SV_Region/paftools_F7_M8/F7
#Find Diploid Assembled Regions: 
grep "^R" LCL7_asm.var.txt|cut -f 2-4 > F7_one_contig.bed
grep "^R" LCL8_asm.var.txt|cut -f 2-4 > M8_one_contig.bed


#Combine Regions:
grep -E -v "chrX|chrY" F7_one_contig.bed > F7_one_contig_autosomes.bed
grep -E -v "chrX|chrY" M8_one_contig.bed > M8_one_contig_autosomes.bed

#Discover Variants: & Compare Variants:
cd /Quartet_data/SV_combine/pacbio_ont_pacbio2_combine_3/Assembly_SV_Region/paftools_F7_M8/F7/Confident_regions/
/software/Jasmine-1.0.1/run.sh  \
max_dist=1000 min_support=0 threads=2 --keep_var_ids --ignore_strand \
file_list=Tier1_paftools.txt out_file=Tier1_paftools.vcf

cd /Quartet_data/SV_combine/pacbio_ont_pacbio2_combine_3/Assembly_SV_Region/paftools_F7_M8/M8/Confident_regions/
/software/Jasmine-1.0.1/run.sh  \
max_dist=1000 min_support=0 threads=2 --keep_var_ids --ignore_strand \
file_list=Tier1_paftools.txt out_file=Tier1_paftools.vcf

#SVwiden
mkdir SVwiden_original
cd SVwiden_original
grep -E "SUPP_VEC=0" ../Tier1_paftools.vcf|grep -E -v "chrX|chrY"|sort -k1,1 -k2,2n > paftools_NoMatchTier1.vcf
/Pacbio_command/vcf2bed.sh paftools_NoMatchTier1.vcf| \
grep -v "#"|sort -k1,1 -k2,2n|cut -f 1-4 -> paftools_NoMatchTier1.bed
ref=/hg38_ref/chr1_22_XYM.fa/chr1_22_XYM.fa
cat <(grep "#" ../Tier1_paftools.vcf) paftools_NoMatchTier1.vcf > paftools_NoMatchTier1_h.vcf

#svanalyzer widen
svanalyzer widen --ref $ref --variants paftools_NoMatchTier1_h.vcf
grep -E -o "REFWIDENED=[a-z0-9:-]+" widened.vcf|wc -l
perl -lane 'if(/REFWIDENED=(.*?):([0-9]+)-([0-9]+)/){print join("\t",$1,$2,$3,$F[2])}' widened.vcf| \
awk -v OFS='\t' '{if($3<=$2){print $1,$3,$2,$4}else{print}}' -| \
bedops --range -50:50 --everything - > widened_50bp.bed

#Trim Regions:
#NoPaftools
bed=/Quartet_data/SV_combine/pacbio_ont_pacbio2_combine_3/Assembly_SV_Region/paftools_F7_M8/F7/F7_one_contig_autosomes.bed
bed=/Quartet_data/SV_combine/pacbio_ont_pacbio2_combine_3/Assembly_SV_Region/paftools_F7_M8/M8/M8_one_contig_autosomes.bed


bedtools subtract -a $bed -b <(cat paftools_NoMatchTier1.bed widened_50bp.bed|sort -k1,1 -k2,2n|bedops --range -50:50 --everything -)| \
sort -k1,1 -k2,2n|bedtools merge -i -|awk '{sum+=$3-$2}END{print "Sum =",sum}'
Sum = 2591967148 #F7
Sum = 2596140552 #M8
bedtools subtract -a $bed -b <(cat paftools_NoMatchTier1.bed widened_50bp.bed|sort -k1,1 -k2,2n|bedops --range -50:50 --everything -)| \
sort -k1,1 -k2,2n|bedtools merge -i -> Tier1_Confident_regions_autosomes_NoPaftools.bed
