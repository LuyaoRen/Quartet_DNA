task rename_fastq {

	File fastq_1_D5
	String fastq_1_D5_name = sub(basename(fastq_1_D5, "_R1.fastq.gz"), "_R1.fq.gz$", "")
	File fastq_1_D6
	String fastq_1_D6_name = sub(basename(fastq_1_D6, "_R1.fastq.gz"), "_R1.fq.gz$", "")
	File fastq_1_F7
	String fastq_1_F7_name = sub(basename(fastq_1_F7, "_R1.fastq.gz"), "_R1.fq.gz$", "")
	File fastq_1_M8
	String fastq_1_M8_name = sub(basename(fastq_1_M8, "_R1.fastq.gz"), "_R1.fq.gz$", "")

	File fastq_2_D5
	String fastq_2_D5_name = sub(basename(fastq_2_D5, "_R2.fastq.gz"), "_R2.fq.gz$", "")
	File fastq_2_D6
	String fastq_2_D6_name = sub(basename(fastq_2_D6, "_R2.fastq.gz"), "_R2.fq.gz$", "")
	File fastq_2_F7
	String fastq_2_F7_name = sub(basename(fastq_2_F7, "_R2.fastq.gz"), "_R2.fq.gz$", "")
	File fastq_2_M8
	String fastq_2_M8_name = sub(basename(fastq_2_M8, "_R2.fastq.gz"), "_R2.fq.gz$", "")

	String project
	String docker
	String cluster_config
	String disk_size

	command <<<
		cp ${fastq_1_D5} ${fastq_1_D5_name}_${project}_LCL5_R1.fastq.gz
		rm ${fastq_1_D5}
		cp ${fastq_1_D6} ${fastq_1_D6_name}_${project}_LCL6_R1.fastq.gz
		rm ${fastq_1_D6}
		cp ${fastq_1_F7} ${fastq_1_F7_name}_${project}_LCL7_R1.fastq.gz
		rm ${fastq_1_F7}
		cp ${fastq_1_M8} ${fastq_1_M8_name}_${project}_LCL8_R1.fastq.gz
		rm ${fastq_1_M8}
		cp ${fastq_2_D5} ${fastq_2_D5_name}_${project}_LCL5_R2.fastq.gz
		rm ${fastq_2_D5}
		cp ${fastq_2_D6} ${fastq_2_D6_name}_${project}_LCL6_R2.fastq.gz
		rm ${fastq_2_D6}
		cp ${fastq_2_F7} ${fastq_2_F7_name}_${project}_LCL7_R2.fastq.gz
		rm ${fastq_2_F7}
		cp ${fastq_2_M8} ${fastq_2_M8_name}_${project}_LCL8_R2.fastq.gz
		rm ${fastq_2_M8}

	>>>

	runtime {
		docker:docker
		cluster:cluster_config
		systemDisk:"cloud_ssd 40"
		dataDisk:"cloud_ssd " + disk_size + " /cromwell_root/"
	}

	output {
		File fastq_1_D5_renamed = "${fastq_1_D5_name}_${project}_LCL5_R1.fastq.gz"
		File fastq_1_D6_renamed = "${fastq_1_D6_name}_${project}_LCL6_R1.fastq.gz"
		File fastq_1_F7_renamed = "${fastq_1_F7_name}_${project}_LCL7_R1.fastq.gz"
		File fastq_1_M8_renamed = "${fastq_1_M8_name}_${project}_LCL8_R1.fastq.gz"
		File fastq_2_D5_renamed = "${fastq_2_D5_name}_${project}_LCL5_R2.fastq.gz"
		File fastq_2_D6_renamed = "${fastq_2_D6_name}_${project}_LCL6_R2.fastq.gz"
		File fastq_2_F7_renamed = "${fastq_2_F7_name}_${project}_LCL7_R2.fastq.gz"
		File fastq_2_M8_renamed = "${fastq_2_M8_name}_${project}_LCL8_R2.fastq.gz"
	}
}