task replicate_consensus {
	File chromo_file
	String chromo
	String docker
	String cluster_config
	String disk_size
	
	command <<<
		cat ${chromo_file} | grep -v '##' > ${chromo}.txt
		python /opt/replicates_consensus.py -vcf ${chromo}.txt -prefix ${chromo}


	>>>

	runtime {
		docker:docker
		cluster: cluster_config
		systemDisk: "cloud_ssd 40"
		dataDisk: "cloud_ssd " + disk_size + " /cromwell_root/"
	}
	output {
		File chromo_consensus = "${chromo}_all_summary.txt"
		File consensus_vcf = "${chromo}_consensus.vcf"
	}
}
