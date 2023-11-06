task extract_tables_vcf {

	File hap
	String project
	String docker
	String cluster_config
	String disk_size

	command <<<
		python /opt/extract_tables.py -hap ${hap} -project ${project}
	>>>

	runtime {
		docker:docker
		cluster:cluster_config
		systemDisk:"cloud_ssd 40"
		dataDisk:"cloud_ssd " + disk_size + " /cromwell_root/"
	}

	output {
		File variant_calling = "variants.calling.qc.txt"
	}
}