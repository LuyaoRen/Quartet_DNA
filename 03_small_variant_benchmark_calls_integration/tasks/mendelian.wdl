task mendelian {
	File consensus_vcf
	File ref_dir
	String fasta
	String chromo
	String docker
	String cluster_config
	String disk_size
	
	command <<<
		export LD_LIBRARY_PATH=/opt/htslib-1.9
		nt=$(nproc)

		echo -e "${chromo}\tLCL8_consensus_call\t0\t0\t2\t-9\n${chromo}\tLCL7_consensus_call\t0\t0\t1\t-9\n${chromo}\tLCL5_consensus_call\tLCL7_consensus_call\tLCL8_consensus_call\t2\t-9" > ${chromo}.D5.ped

		mkdir VBT_D5
		/opt/VBT-TrioAnalysis/vbt mendelian -ref ${ref_dir}/${fasta} -mother ${consensus_vcf} -father ${consensus_vcf} -child ${consensus_vcf} -pedigree ${chromo}.D5.ped -outDir VBT_D5 -out-prefix ${chromo}.D5 --output-violation-regions -thread-count $nt

		cat VBT_D5/${chromo}.D5_trio.vcf > ${chromo}.D5.vcf

		echo -e "${chromo}\tLCL8_consensus_call\t0\t0\t2\t-9\n${chromo}\tLCL7_consensus_call\t0\t0\t1\t-9\n${chromo}\tLCL6_consensus_call\tLCL7_consensus_call\tLCL8_consensus_call\t2\t-9" > ${chromo}.D6.ped

		mkdir VBT_D6
		/opt/VBT-TrioAnalysis/vbt mendelian -ref ${ref_dir}/${fasta} -mother ${consensus_vcf} -father ${consensus_vcf} -child ${consensus_vcf} -pedigree ${chromo}.D6.ped -outDir VBT_D6 -out-prefix ${chromo}.D6 --output-violation-regions -thread-count $nt

		cat VBT_D6/${chromo}.D6_trio.vcf > ${chromo}.D6.vcf
	>>>

	runtime {
		docker:docker
		cluster: cluster_config
		systemDisk: "cloud_ssd 40"
		dataDisk: "cloud_ssd " + disk_size + " /cromwell_root/"
	}
	output {
		File D5_ped = "${chromo}.D5.ped"
		File D6_ped = "${chromo}.D6.ped"
		Array[File] D5_mendelian = glob("VBT_D5/*")
		Array[File] D6_mendelian = glob("VBT_D6/*")
		File D5_trio_vcf = "${chromo}.D5.vcf"
		File D6_trio_vcf = "${chromo}.D6.vcf"
		File family_vcf = "${chromo}.vcf"
	}
}



