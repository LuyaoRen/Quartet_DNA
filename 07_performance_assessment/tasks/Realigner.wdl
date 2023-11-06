task Realigner {

    File ref_dir
    File dbmills_dir

	String SENTIEON_INSTALL_DIR
	String fasta

	File Dedup_bam
	File Dedup_bam_index
	String sample = basename(Dedup_bam,".sorted.deduped.bam")
	String db_mills
	String docker	
	String cluster_config
	String disk_size


	command <<<
	set -o pipefail
	set -e
	export SENTIEON_LICENSE=192.168.0.55:8990
	nt=$(nproc)
	${SENTIEON_INSTALL_DIR}/bin/sentieon driver -r ${ref_dir}/${fasta} -t $nt -i ${Dedup_bam} --algo Realigner -k ${dbmills_dir}/${db_mills} ${sample}.sorted.deduped.realigned.bam
	
	>>>

	runtime {
		docker:docker
		cluster: cluster_config
    	systemDisk: "cloud_ssd 40"
    	dataDisk: "cloud_ssd " + disk_size + " /cromwell_root/"
	}

	output {
		File realigner_bam = "${sample}.sorted.deduped.realigned.bam"
		File realigner_bam_index = "${sample}.sorted.deduped.realigned.bam.bai"

	}
}


