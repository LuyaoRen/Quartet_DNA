task mapping {

    File ref_dir
    String fasta
	File fastq_1
	File fastq_2

	String SENTIEON_INSTALL_DIR
	String SENTIEON_LICENSE
	String group
	String sample
	String project
	String pl
	String user_define_name = sub(basename(fastq_1, "_R1.fastq.gz"), "_R1.fq.gz$", "")
	String docker
	String cluster_config
	String disk_size

	command <<<
		set -o pipefail
		set -e	
		export SENTIEON_LICENSE=${SENTIEON_LICENSE}
		nt=$(nproc)
		${SENTIEON_INSTALL_DIR}/bin/bwa mem -M -R "@RG\tID:${group}\tSM:${sample}\tPL:${pl}" -t $nt -K 10000000 ${ref_dir}/${fasta} ${fastq_1} ${fastq_2} | ${SENTIEON_INSTALL_DIR}/bin/sentieon util sort -o ${user_define_name}_${project}_${sample}.sorted.bam -t $nt --sam2bam -i -
	>>>

	runtime {
		docker:docker
    	cluster: cluster_config
    	systemDisk: "cloud_ssd 40"
    	dataDisk: "cloud_ssd " + disk_size + " /cromwell_root/"
	}
	output {
		File sorted_bam = "${user_define_name}_${project}_${sample}.sorted.bam"
		File sorted_bam_index = "${user_define_name}_${project}_${sample}.sorted.bam.bai"
	}
}
