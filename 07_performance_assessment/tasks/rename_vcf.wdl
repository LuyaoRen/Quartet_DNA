task rename_vcf {

	File vcf_D5
	String vcf_D5_name = basename(vcf_D5, ".vcf")
	File vcf_D6
	String vcf_D6_name = basename(vcf_D6, ".vcf")
	File vcf_F7
	String vcf_F7_name = basename(vcf_F7, ".vcf")
	File vcf_M8
	String vcf_M8_name = basename(vcf_M8, ".vcf")

	String project
	String docker
	String cluster_config
	String disk_size

	command <<<
		cp ${vcf_D5} ${vcf_D5_name}_${project}_LCL5.vcf
		rm ${vcf_D5}
		cp ${vcf_D6} ${vcf_D6_name}_${project}_LCL6.vcf
		rm ${vcf_D6}
		cp ${vcf_F7} ${vcf_F7_name}_${project}_LCL7.vcf
		rm ${vcf_F7}
		cp ${vcf_M8} ${vcf_M8_name}_${project}_LCL8.vcf
		rm ${vcf_M8}
	>>>

	runtime {
		docker:docker
		cluster:cluster_config
		systemDisk:"cloud_ssd 40"
		dataDisk:"cloud_ssd " + disk_size + " /cromwell_root/"
	}

	output {
		File vcf_D5_renamed = "${vcf_D5_name}_${project}_LCL5.vcf"
		File vcf_D6_renamed = "${vcf_D6_name}_${project}_LCL6.vcf"
		File vcf_F7_renamed = "${vcf_F7_name}_${project}_LCL7.vcf"
		File vcf_M8_renamed = "${vcf_M8_name}_${project}_LCL8.vcf"
	}
}