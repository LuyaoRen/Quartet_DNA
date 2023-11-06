task CallableLoci {
	
	File bam
	File bam_index
	File ref_dir
	String fasta
	String sample
	String docker
	String disk_size
	String cluster_config
	

	command <<<
		set -o pipefail
		set -e
		java -jar /usr/GenomeAnalysisTK.jar \
		-T CallableLoci \
		-R ${ref_dir}/${fasta} \
		-I ${bam} \
		--maxDepth 300 \
		--maxFractionOfReadsWithLowMAPQ 0.1 \
		--maxLowMAPQ 1 \
		--minBaseQuality 20 \
		--minMappingQuality 20 \
		--minDepth 10 \
		--minDepthForLowMAPQ 10 \
		-summary ${sample}_table.txt \
		-o ${sample}_callable_status.bed
		>>>

		runtime {
		docker:docker
    	cluster:cluster_config
    	systemDisk: "cloud_ssd 40"
    	dataDisk: "cloud_ssd " + disk_size + " /cromwell_root/"
	}

	output {
		File summary = "${sample}_table.txt"
		File bed = "${sample}_callable_status.bed"
	}
}




