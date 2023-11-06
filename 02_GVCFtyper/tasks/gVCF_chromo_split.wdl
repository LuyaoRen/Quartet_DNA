task gVCF_chromo_split {
	
	File gvcf
	File gvcf_idx
	String SENTIEON_INSTALL_DIR
	String sample
	String docker
	String cluster_config
	String disk_size
	

	command <<<
		cat ${gvcf} | grep '#' > header
		for i in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX
		do
			cat ${gvcf} | grep -v '#' | grep -w $i | cat header -  > ${sample}.$i.vcf
		done

		set -o pipefail
		set -e
		export SENTIEON_LICENSE=192.168.0.55:8990
		nt=$(nproc)	
		${SENTIEON_INSTALL_DIR}/bin/sentieon util vcfindex ${sample}.chr1.vcf
		${SENTIEON_INSTALL_DIR}/bin/sentieon util vcfindex ${sample}.chr2.vcf
		${SENTIEON_INSTALL_DIR}/bin/sentieon util vcfindex ${sample}.chr3.vcf
		${SENTIEON_INSTALL_DIR}/bin/sentieon util vcfindex ${sample}.chr4.vcf
		${SENTIEON_INSTALL_DIR}/bin/sentieon util vcfindex ${sample}.chr5.vcf
		${SENTIEON_INSTALL_DIR}/bin/sentieon util vcfindex ${sample}.chr6.vcf
		${SENTIEON_INSTALL_DIR}/bin/sentieon util vcfindex ${sample}.chr7.vcf
		${SENTIEON_INSTALL_DIR}/bin/sentieon util vcfindex ${sample}.chr8.vcf
		${SENTIEON_INSTALL_DIR}/bin/sentieon util vcfindex ${sample}.chr9.vcf
		${SENTIEON_INSTALL_DIR}/bin/sentieon util vcfindex ${sample}.chr10.vcf
		${SENTIEON_INSTALL_DIR}/bin/sentieon util vcfindex ${sample}.chr11.vcf
		${SENTIEON_INSTALL_DIR}/bin/sentieon util vcfindex ${sample}.chr12.vcf
		${SENTIEON_INSTALL_DIR}/bin/sentieon util vcfindex ${sample}.chr13.vcf
		${SENTIEON_INSTALL_DIR}/bin/sentieon util vcfindex ${sample}.chr14.vcf
		${SENTIEON_INSTALL_DIR}/bin/sentieon util vcfindex ${sample}.chr15.vcf
		${SENTIEON_INSTALL_DIR}/bin/sentieon util vcfindex ${sample}.chr16.vcf
		${SENTIEON_INSTALL_DIR}/bin/sentieon util vcfindex ${sample}.chr17.vcf
		${SENTIEON_INSTALL_DIR}/bin/sentieon util vcfindex ${sample}.chr18.vcf
		${SENTIEON_INSTALL_DIR}/bin/sentieon util vcfindex ${sample}.chr19.vcf
		${SENTIEON_INSTALL_DIR}/bin/sentieon util vcfindex ${sample}.chr20.vcf
		${SENTIEON_INSTALL_DIR}/bin/sentieon util vcfindex ${sample}.chr21.vcf
		${SENTIEON_INSTALL_DIR}/bin/sentieon util vcfindex ${sample}.chr22.vcf
		${SENTIEON_INSTALL_DIR}/bin/sentieon util vcfindex ${sample}.chrX.vcf


		>>>

		runtime {
		docker:docker
    	cluster:cluster_config
    	systemDisk: "cloud_ssd 40"
    	dataDisk: "cloud_ssd " + disk_size + " /cromwell_root/"
	}

	output {
		File chr1_gvcf = "${sample}.chr1.vcf"
		File chr2_gvcf = "${sample}.chr2.vcf"
		File chr3_gvcf = "${sample}.chr3.vcf"
		File chr4_gvcf = "${sample}.chr4.vcf"
		File chr5_gvcf = "${sample}.chr5.vcf"
		File chr6_gvcf = "${sample}.chr6.vcf"
		File chr7_gvcf = "${sample}.chr7.vcf"
		File chr8_gvcf = "${sample}.chr8.vcf"
		File chr9_gvcf = "${sample}.chr9.vcf"
		File chr10_gvcf = "${sample}.chr10.vcf"
		File chr11_gvcf = "${sample}.chr11.vcf"
		File chr12_gvcf = "${sample}.chr12.vcf"
		File chr13_gvcf = "${sample}.chr13.vcf"
		File chr14_gvcf = "${sample}.chr14.vcf"
		File chr15_gvcf = "${sample}.chr15.vcf"
		File chr16_gvcf = "${sample}.chr16.vcf"
		File chr17_gvcf = "${sample}.chr17.vcf"
		File chr18_gvcf = "${sample}.chr18.vcf"
		File chr19_gvcf = "${sample}.chr19.vcf"
		File chr20_gvcf = "${sample}.chr20.vcf"
		File chr21_gvcf = "${sample}.chr21.vcf"
		File chr22_gvcf = "${sample}.chr22.vcf"
		File chrX_gvcf = "${sample}.chrX.vcf"
		File chr1_gvcf_idx = "${sample}.chr1.vcf.idx"
		File chr2_gvcf_idx = "${sample}.chr2.vcf.idx"
		File chr3_gvcf_idx = "${sample}.chr3.vcf.idx"
		File chr4_gvcf_idx = "${sample}.chr4.vcf.idx"
		File chr5_gvcf_idx = "${sample}.chr5.vcf.idx"
		File chr6_gvcf_idx = "${sample}.chr6.vcf.idx"
		File chr7_gvcf_idx = "${sample}.chr7.vcf.idx"
		File chr8_gvcf_idx = "${sample}.chr8.vcf.idx"
		File chr9_gvcf_idx = "${sample}.chr9.vcf.idx"
		File chr10_gvcf_idx = "${sample}.chr10.vcf.idx"
		File chr11_gvcf_idx = "${sample}.chr11.vcf.idx"
		File chr12_gvcf_idx = "${sample}.chr12.vcf.idx"
		File chr13_gvcf_idx = "${sample}.chr13.vcf.idx"
		File chr14_gvcf_idx = "${sample}.chr14.vcf.idx"
		File chr15_gvcf_idx = "${sample}.chr15.vcf.idx"
		File chr16_gvcf_idx = "${sample}.chr16.vcf.idx"
		File chr17_gvcf_idx = "${sample}.chr17.vcf.idx"
		File chr18_gvcf_idx = "${sample}.chr18.vcf.idx"
		File chr19_gvcf_idx = "${sample}.chr19.vcf.idx"
		File chr20_gvcf_idx = "${sample}.chr20.vcf.idx"
		File chr21_gvcf_idx = "${sample}.chr21.vcf.idx"
		File chr22_gvcf_idx = "${sample}.chr22.vcf.idx"
		File chrX_gvcf_idx = "${sample}.chrX.vcf.idx"
	}
}




