task joint_variant_calling {
	
	File merged_bed
	String sample
	String docker
	String disk_size
	String cluster_config
	

	command <<<
		
		sentieon driver -r REFERENCE --algo GVCFtyper -v s1_VARIANT_GVCF \ -v s2_VARIANT_GVCF -v s3_VARIANT_GVCF VARIANT_VCF

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




