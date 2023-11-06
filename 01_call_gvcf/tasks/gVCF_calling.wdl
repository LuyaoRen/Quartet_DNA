task gVCF_calling {
	
	File bam
	File bai
	File ref_dir
	String SENTIEON_INSTALL_DIR
	String fasta
	String sample
	String docker
	String disk_size
	String cluster_config
	

	command <<<
		export SENTIEON_LICENSE=192.168.0.55:8990
		nt=$(nproc)
		${SENTIEON_INSTALL_DIR}/bin/sentieon driver -r ${ref_dir}/${fasta} -t $nt -i ${bam} -q ${sample}_RECAL_DATA_TABLE --algo Haplotyper --emit_mode gvcf \ ${sample}_VARIANT_GVCF

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




