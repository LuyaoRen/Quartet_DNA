task CallableLoci {
	
	File bed
	String sample
	String docker
	String disk_size
	String cluster_config
	

	command <<<

		cat ${bed} | grep CALLABLE > ${sample}.CALLABLE.bed
		>>>

		runtime {
		docker:docker
    	cluster:cluster_config
    	systemDisk: "cloud_ssd 40"
    	dataDisk: "cloud_ssd " + disk_size + " /cromwell_root/"
	}

	output {
		File callable_bed = "${sample}.CALLABLE.bed"
	}
}




