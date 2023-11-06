task merge_sentieon_metrics {
	Array[File] quality_yield_header
	Array[File] wgs_metrics_algo_header
	Array[File] aln_metrics_header
	Array[File] is_metrics_header

	Array[File] quality_yield_data
	Array[File] wgs_metrics_algo_data
	Array[File] aln_metrics_data
	Array[File] is_metrics_data

	String project
	String docker
	String cluster_config
	String disk_size
	
	command <<<
		echo '''Sample''' > sample_column

		cat ${sep=" " quality_yield_header} | sed -n '1,1p' | cat - ${sep=" " quality_yield_data} > quality_yield_all
		ls ${sep=" " quality_yield_data} | cut -d '.' -f1 | cat sample_column - | paste - quality_yield_all > ${project}.quality_yield.txt

		cat ${sep=" " wgs_metrics_algo_header} | sed -n '1,1p' | cat - ${sep=" " wgs_metrics_algo_data} > wgs_metrics_all	
		ls ${sep=" " wgs_metrics_algo_data} | cut -d '.' -f1 | cat sample_column - | paste - wgs_metrics_all > ${project}.wgs_metrics_data.txt
		
		cat ${sep=" " aln_metrics_header} | sed -n '1,1p' | cat - ${sep=" " aln_metrics_data} > aln_metrics_all
		ls ${sep=" " aln_metrics_data} | cut -d '.' -f1 | cat sample_column - | paste - aln_metrics_all > ${project}.aln_metrics.txt

		cat ${sep=" " is_metrics_header} | sed -n '1,1p' | cat - ${sep=" " is_metrics_data} > is_metrics_all
		ls ${sep=" " is_metrics_data} | cut -d '.' -f1 | cat sample_column - | paste - is_metrics_all > ${project}.is_metrics.txt

	>>>

	runtime {
		docker:docker
		cluster: cluster_config
		systemDisk: "cloud_ssd 40"
		dataDisk: "cloud_ssd " + disk_size + " /cromwell_root/"
	}

	output {
		File quality_yield_summary = "${project}.quality_yield.txt"
		File wgs_metrics_summary = "${project}.wgs_metrics_data.txt"
		File aln_metrics_summary = "${project}.aln_metrics.txt"
		File is_metrics_summary = "${project}.is_metrics.txt"
	}
}