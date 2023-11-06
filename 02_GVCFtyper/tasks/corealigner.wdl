task corealigner {
	
    File ref_dir
    File dbsnp_dir
    File dbmills_dir

	String sample
	String SENTIEON_INSTALL_DIR
	String docker
	String cluster_config
	String fasta

	String dbsnp
	String db_mills
	File tumor_recaled_bam
	File tumor_recaled_bam_index
	File normal_recaled_bam
	File normal_recaled_bam_index
	String disk_size


	command <<<
		set -o pipefail
		set -e
		export SENTIEON_LICENSE=192.168.0.55:8990
		nt=$(nproc)
		${SENTIEON_INSTALL_DIR}/bin/sentieon driver -r ${ref_dir}/${fasta} -t $nt -i ${tumor_recaled_bam} -i ${normal_recaled_bam} --algo Realigner -k ${db_mills} -k ${dbsnp} ${sample}_corealigned.bam				
	>>>
	
	runtime {
		docker:docker
    	cluster: cluster_config
    	systemDisk: "cloud_ssd 40"
    	dataDisk: "cloud_ssd " + disk_size + " /cromwell_root/"
	}

	output {
		File corealigner_bam = "${sample}_corealigned.bam"
		File corealigner_bam_index = "${sample}_corealigned.bam.bai"
	}	
}



