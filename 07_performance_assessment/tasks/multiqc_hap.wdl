task multiqc_hap {

	Array[File] summary

	String docker
	String cluster_config
	String disk_size

	command <<<
		set -o pipefail
		set -e
		mkdir -p /cromwell_root/tmp/benchmark
		cp ${sep=" " summary} /cromwell_root/tmp/benchmark
		multiqc /cromwell_root/tmp/
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
		File hap = "multiqc_happy_data.json"
	}
}