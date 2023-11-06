task mergeBed {
	
	Array[File] callable_bed
	String sample
	String docker
	String disk_size
	String cluster_config
	

	command <<<
		
		/opt/ccdg/bedtools-2.27.1/bin/bedtools multiinter -i ${sep=" " callable_bed} > ${sample}.CALLABLE.merged.bed

		>>>

		runtime {
		docker:docker
    	cluster:cluster_config
    	systemDisk: "cloud_ssd 40"
    	dataDisk: "cloud_ssd " + disk_size + " /cromwell_root/"
	}

	output {
		File merged_bed = "${sample}.CALLABLE.merged.bed"
	}
}




