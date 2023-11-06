task Haplotyper {
	
	
    File ref_dir
	String SENTIEON_INSTALL_DIR
	String fasta
	File recaled_bam
	File recaled_bam_index
	String sample
	String docker
	String cluster_config
	String disk_size

command <<<
		set -o pipefail
		set -e
		export SENTIEON_LICENSE=192.168.0.55:8990
		nt=$(nproc)	
		${SENTIEON_INSTALL_DIR}/bin/sentieon driver -r ${ref_dir}/${fasta} -t $nt -i ${recaled_bam} --algo Haplotyper --emit_mode gvcf ${sample}_hc.g.vcf
	>>>
	
	runtime {
		docker:docker
    	cluster: cluster_config
    	systemDisk: "cloud_ssd 40"
    	dataDisk: "cloud_ssd " + disk_size + " /cromwell_root/"
	}

	output {
		File vcf = "${sample}_hc.g.vcf"
		File vcf_idx = "${sample}_hc.g.vcf.idx"
	}
}


