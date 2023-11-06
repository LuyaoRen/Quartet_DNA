task multiqc {

	Array[File] read1_zip
	Array[File] read2_zip

	Array[File] txt1
	Array[File] txt2

	Array[File] summary

	String docker
	String cluster_config
	String disk_size

	command <<<
		set -o pipefail
		set -e
		mkdir -p /cromwell_root/tmp/fastqc
		mkdir -p /cromwell_root/tmp/fastqscreen
		mkdir -p /cromwell_root/tmp/benchmark

		cp ${sep=" " read1_zip} ${sep=" " read2_zip} /cromwell_root/tmp/fastqc
		cp ${sep=" " txt1} ${sep=" " txt2} /cromwell_root/tmp/fastqscreen
		cp ${sep=" " summary} /cromwell_root/tmp/benchmark

		multiqc /cromwell_root/tmp/

		cat multiqc_data/multiqc_general_stats.txt > multiqc_general_stats.txt
		cat multiqc_data/multiqc_fastq_screen.txt > multiqc_fastq_screen.txt
		cat multiqc_data/multiqc_happy_data.json > multiqc_happy_data.json
	
	>>>

	runtime {
		docker:docker
		cluster:cluster_config
		systemDisk:"cloud_ssd 40"
		dataDisk:"cloud_ssd " + disk_size + " /cromwell_root/"
	}

	output {
		File multiqc_html = "multiqc_report.html"
		Array[File] multiqc_txt = glob("multiqc_data/*")
		File? fastqc = "multiqc_general_stats.txt"
		File? fastqscreen = "multiqc_fastq_screen.txt"
		File hap = "multiqc_happy_data.json"
	}
}