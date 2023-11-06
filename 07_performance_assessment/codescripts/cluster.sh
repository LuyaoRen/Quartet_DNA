cat benchmark.men.vote.diffbed.lengthlessthan50.txt | awk '{print $1"\t"$2"\t"".""\t"$35"\t"$7"\t.\t.\t.\tGT\t"$6}' | grep -v '0/0' > LCL5.body
cat benchmark.men.vote.diffbed.lengthlessthan50.txt | awk '{print $1"\t"$2"\t"".""\t"$35"\t"$15"\t.\t.\t.\tGT\t"$14}' | grep -v '0/0' > LCL6.body
cat benchmark.men.vote.diffbed.lengthlessthan50.txt | awk '{print $1"\t"$2"\t"".""\t"$35"\t"$23"\t.\t.\t.\tGT\t"$22}' | grep -v '0/0'> LCL7.body
cat benchmark.men.vote.diffbed.lengthlessthan50.txt | awk '{print $1"\t"$2"\t"".""\t"$35"\t"$31"\t.\t.\t.\tGT\t"$30}'| grep -v '0/0' > LCL8.body
cat header5 LCL5.body > LCL5.beforediffbed.vcf
cat header6 LCL6.body > LCL6.beforediffbed.vcf
cat header7 LCL7.body > LCL7.beforediffbed.vcf
cat header8 LCL8.body > LCL8.beforediffbed.vcf
rtg bgzip *beforediffbed.vcf
rtg index *beforediffbed.vcf.gz

rtg vcffilter -i LCL5.beforediffbed.vcf.gz --exclude-bed=/mnt/pgx_src_data_pool_4/home/renluyao/manuscript/benchmark_calls/MIE/diff.merged.bed -o LCL5.afterfilterdiffbed.vcf.gz
rtg vcffilter -i LCL6.beforediffbed.vcf.gz --exclude-bed=/mnt/pgx_src_data_pool_4/home/renluyao/manuscript/benchmark_calls/MIE/diff.merged.bed -o LCL6.afterfilterdiffbed.vcf.gz
rtg vcffilter -i LCL7.beforediffbed.vcf.gz --exclude-bed=/mnt/pgx_src_data_pool_4/home/renluyao/manuscript/benchmark_calls/MIE/diff.merged.bed -o LCL7.afterfilterdiffbed.vcf.gz
rtg vcffilter -i LCL8.beforediffbed.vcf.gz --exclude-bed=/mnt/pgx_src_data_pool_4/home/renluyao/manuscript/benchmark_calls/MIE/diff.merged.bed -o LCL8.afterfilterdiffbed.vcf.gz

/mnt/pgx_src_data_pool_4/home/renluyao/softwares/annovar/table_annovar.pl LCL5.beforediffbed.vcf.gz /mnt/pgx_src_data_pool_4/home/renluyao/softwares/annovar/humandb \
-buildver hg38 \
-out LCL5 \
-remove \
-protocol 1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,clinvar_20190305,gnomad211_genome \
-operation f,f,f,f,f,f,f,f \
-nastring . \
-vcfinput \
--thread 8

rtg vcfeval -b /mnt/pgx_src_data_pool_4/home/renluyao/Quartet/GIAB/NA12878_HG001/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz -c LCL5.afterfilterdiffbed.vcf.gz -o LCL5_NIST -t /mnt/pgx_src_data_pool_4/home/renluyao/annotation/hg38/GRCh38.d1.vd1.sdf/
rtg vcfeval -b /mnt/pgx_src_data_pool_4/home/renluyao/Quartet/GIAB/NA12878_HG001/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz -c LCL6.afterfilterdiffbed.vcf.gz -o LCL6_NIST -t /mnt/pgx_src_data_pool_4/home/renluyao/annotation/hg38/GRCh38.d1.vd1.sdf/
rtg vcfeval -b /mnt/pgx_src_data_pool_4/home/renluyao/Quartet/GIAB/NA12878_HG001/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz -c LCL7.afterfilterdiffbed.vcf.gz -o LCL7_NIST -t /mnt/pgx_src_data_pool_4/home/renluyao/annotation/hg38/GRCh38.d1.vd1.sdf/
rtg vcfeval -b /mnt/pgx_src_data_pool_4/home/renluyao/Quartet/GIAB/NA12878_HG001/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz -c LCL8.afterfilterdiffbed.vcf.gz -o LCL8_NIST -t /mnt/pgx_src_data_pool_4/home/renluyao/annotation/hg38/GRCh38.d1.vd1.sdf/


zcat LCL5.afterfilterdiffbed.vcf.gz | grep -v '#' | awk '{ if ((length($4) == 1) && (length($5) == 1)) { print } }' | wc -l
zcat LCL6.afterfilterdiffbed.vcf.gz | grep -v '#' | awk '{ if ((length($4) == 1) && (length($5) == 1)) { print } }' | wc -l
zcat LCL7.afterfilterdiffbed.vcf.gz | grep -v '#' | awk '{ if ((length($4) == 1) && (length($5) == 1)) { print } }' | wc -l
zcat LCL8.afterfilterdiffbed.vcf.gz | grep -v '#' | awk '{ if ((length($4) == 1) && (length($5) == 1)) { print } }' | wc -l


zcat LCL5.afterfilterdiffbed.vcf.gz | grep -v '#' | awk '{ if ((length($4) == 1) && (length($5) < 11) && (length($5) > 1)) { print } }' | wc -l
zcat LCL6.afterfilterdiffbed.vcf.gz | grep -v '#' | awk '{ if ((length($4) == 1) && (length($5) < 11) && (length($5) > 1)) { print } }' | wc -l
zcat LCL7.afterfilterdiffbed.vcf.gz | grep -v '#' | awk '{ if ((length($4) == 1) && (length($5) < 11) && (length($5) > 1)) { print } }' | wc -l
zcat LCL8.afterfilterdiffbed.vcf.gz | grep -v '#' | awk '{ if ((length($4) == 1) && (length($5) < 11) && (length($5) > 1)) { print } }' | wc -l


zcat LCL5.afterfilterdiffbed.vcf.gz | grep -v '#' | awk '{ if ((length($4) == 1) && (length($5) > 10)) { print } }' | wc -l
zcat LCL6.afterfilterdiffbed.vcf.gz | grep -v '#' | awk '{ if ((length($4) == 1) && (length($5) > 10)) { print } }' | wc -l
zcat LCL7.afterfilterdiffbed.vcf.gz | grep -v '#' | awk '{ if ((length($4) == 1) && (length($5) > 10)) { print } }' | wc -l
zcat LCL8.afterfilterdiffbed.vcf.gz | grep -v '#' | awk '{ if ((length($4) == 1) && (length($5) > 10)) { print } }' | wc -l

bedtools subtract -a LCL5.27.homo_ref.consensus.bed -b /mnt/pgx_src_data_pool_4/home/renluyao/manuscript/benchmark_calls/MIE/diff.merged.bed > LCL5.27.homo_ref.consensus.filtereddiffbed.bed
bedtools subtract -a LCL6.27.homo_ref.consensus.bed -b /mnt/pgx_src_data_pool_4/home/renluyao/manuscript/benchmark_calls/MIE/diff.merged.bed > LCL6.27.homo_ref.consensus.filtereddiffbed.bed
bedtools subtract -a LCL7.27.homo_ref.consensus.bed -b /mnt/pgx_src_data_pool_4/home/renluyao/manuscript/benchmark_calls/MIE/diff.merged.bed > LCL7.27.homo_ref.consensus.filtereddiffbed.bed
bedtools subtract -a LCL8.27.homo_ref.consensus.bed -b /mnt/pgx_src_data_pool_4/home/renluyao/manuscript/benchmark_calls/MIE/diff.merged.bed > LCL8.27.homo_ref.consensus.filtereddiffbed.bed


python vcf2bed.py LCL5.body LCL5.variants.bed
python vcf2bed.py LCL6.body LCL6.variants.bed
python vcf2bed.py LCL7.body LCL7.variants.bed
python vcf2bed.py LCL8.body LCL8.variants.bed



cat /mnt/pgx_src_data_pool_4/home/renluyao/manuscript/benchmark_calls/all_info/LCL5.variants.bed | cut -f1,11,12 | cat - LCL5.27.homo_ref.consensus.filtereddiffbed.bed | sort -k1,1 -k2,2n > LCL5.high.confidence.bed




