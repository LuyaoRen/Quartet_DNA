task GVCFtyper {
	
    File ref_dir
	String SENTIEON_INSTALL_DIR
	String fasta
	Array[File] vcf
	Array[File] vcf_idx
	String quartet_sample
	String docker
	String cluster_config
	String disk_size

command <<<
		set -o pipefail
		set -e
		export SENTIEON_LICENSE=192.168.0.55:8990
		nt=$(nproc)	
		${SENTIEON_INSTALL_DIR}/bin/sentieon driver -r ${ref_dir}/${fasta} --algo GVCFtyper ${quartet_sample}.joint.g.vcf ${sep=" " vcf}
	>>>
	
	runtime {
		docker:docker
    	cluster: cluster_config
    	systemDisk: "cloud_ssd 40"
    	dataDisk: "cloud_ssd " + disk_size + " /cromwell_root/"
	}

	output {
		File gvcf = "${quartet_sample}.joint.g.vcf"
		File gvcf_idx = "${quartet_sample}.joint.g.vcf.idx"
	}
}


