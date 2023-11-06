task Dedup {

	String SENTIEON_INSTALL_DIR
	File sorted_bam
	File sorted_bam_index
	String sample = basename(sorted_bam,".sorted.bam")
	String docker
	String cluster_config
	String disk_size


	command <<<
		set -o pipefail
		set -e
		export SENTIEON_LICENSE=192.168.0.55:8990
		nt=$(nproc)
		${SENTIEON_INSTALL_DIR}/bin/sentieon driver -t $nt -i ${sorted_bam} --algo LocusCollector --fun score_info ${sample}_score.txt
		${SENTIEON_INSTALL_DIR}/bin/sentieon driver -t $nt -i ${sorted_bam} --algo Dedup --rmdup --score_info ${sample}_score.txt --metrics ${sample}_dedup_metrics.txt ${sample}.sorted.deduped.bam	
	>>>
	runtime {
		docker:docker
    	cluster: cluster_config
    	systemDisk: "cloud_ssd 40"
    	dataDisk: "cloud_ssd " + disk_size + " /cromwell_root/"
	}

	output {
		File score = "${sample}_score.txt"
		File dedup_metrics = "${sample}_dedup_metrics.txt"
		File Dedup_bam = "${sample}.sorted.deduped.bam"
		File Dedup_bam_index = "${sample}.sorted.deduped.bam.bai"
	}
}






