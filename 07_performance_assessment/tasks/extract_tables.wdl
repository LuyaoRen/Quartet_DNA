task extract_tables {

	File quality_yield_summary
	File wgs_metrics_summary
	File aln_metrics_summary
	File is_metrics_summary
	File hap
	File fastqc
	File fastqscreen


	String project
	String docker
	String cluster_config
	String disk_size

	command <<<
		python /opt/extract_tables.py -quality ${quality_yield_summary} -depth ${wgs_metrics_summary} -aln ${aln_metrics_summary} -is ${is_metrics_summary} -fastqc ${fastqc} -fastqscreen ${fastqscreen} -hap ${hap} -project ${project}
	>>>

	runtime {
		docker:docker
		cluster:cluster_config
		systemDisk:"cloud_ssd 40"
		dataDisk:"cloud_ssd " + disk_size + " /cromwell_root/"
	}

	output {
		File pre_alignment = "pre_alignment.txt"
		File post_alignment = "post_alignment.txt"
		File variant_calling = "variants.calling.qc.txt"
	}
}