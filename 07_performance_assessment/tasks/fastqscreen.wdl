task fastq_screen {
	File read1
	File read2
	File screen_ref_dir
	File fastq_screen_conf
	String docker
	String project
	String user_define_name = sub(basename(read1, "_R1.fastq.gz"), "_R1.fq.gz$", "")
	String sample
	String cluster_config
	String disk_size

	command <<<
		set -o pipefail
		set -e
		nt=$(nproc)
		mkdir -p /cromwell_root/tmp
		cp -r ${screen_ref_dir} /cromwell_root/tmp/
		cp ${read1} ${user_define_name}_${project}_${sample}_R1.fastq.gz
		cp ${read2} ${user_define_name}_${project}_${sample}_R2.fastq.gz
		fastq_screen --aligner bowtie2 --conf ${fastq_screen_conf} --subset 1000000 --threads $nt ${user_define_name}_${project}_${sample}_R1.fastq.gz
		fastq_screen --aligner bowtie2 --conf ${fastq_screen_conf} --subset 1000000 --threads $nt ${user_define_name}_${project}_${sample}_R2.fastq.gz
	>>>

	runtime {
		docker:docker
		cluster: cluster_config
		systemDisk: "cloud_ssd 40"
		dataDisk: "cloud_ssd " + disk_size + " /cromwell_root/"
	}
	
	output {
		File png1 = "${user_define_name}_${project}_${sample}_R1_screen.png"
		File txt1 = "${user_define_name}_${project}_${sample}_R1_screen.txt"
		File html1 = "${user_define_name}_${project}_${sample}_R1_screen.html"
		File png2 = "${user_define_name}_${project}_${sample}_R2_screen.png"
		File txt2 = "${user_define_name}_${project}_${sample}_R2_screen.txt"
		File html2 = "${user_define_name}_${project}_${sample}_R2_screen.html"
	}
}