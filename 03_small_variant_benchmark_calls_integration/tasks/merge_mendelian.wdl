task merge_mendelian {
	File D5_trio_vcf
	File D6_trio_vcf
	File consensus_vcf
	String chromo
	String docker
	String cluster_config
	String disk_size
	
	command <<<
		cat ${D5_trio_vcf} | grep -v '##' > ${chromo}.D5.txt
		cat ${D6_trio_vcf} | grep -v '##' > ${chromo}.D6.txt
		cat ${consensus_vcf} | grep -v '##' > ${chromo}.consensus.txt
		python /opt/merge_two_family_with_genotype.py -LCL5 ${chromo}.D5.txt -LCL6 ${chromo}.D6.txt -genotype ${chromo}.consensus.txt -family ${chromo}.mendelian
	>>>

	runtime {
		docker:docker
		cluster: cluster_config
		systemDisk: "cloud_ssd 40"
		dataDisk: "cloud_ssd " + disk_size + " /cromwell_root/"
	}
	output {
		File chromo_mendelian = "${chromo}.mendelian.txt"
	}
}

## change mendelian.txt to vcf
