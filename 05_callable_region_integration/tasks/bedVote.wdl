task bedVote {
	
	File merged_bed
	String sample
	String docker
	String disk_size
	String cluster_config
	

	command <<<
		
		python /opt/callable_bed_voting.py -bed ${merged_bed} -prefix ${sample}

		>>>

		runtime {
		docker:docker
    	cluster:cluster_config
    	systemDisk: "cloud_ssd 40"
    	dataDisk: "cloud_ssd " + disk_size + " /cromwell_root/"
	}

	output {
		File consensus_bed = "${sample}.27consensus.bed"
		File filtered_bed = "${sample}.filtered.bed"
	}
}




