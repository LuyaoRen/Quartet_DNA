task sentieon {
	File quality_yield
	File wgs_metrics_algo
	File aln_metrics
	File is_metrics

	String sample = basename(quality_yield,"_deduped_QualityYield.txt")
	String docker
	String cluster_config
	String disk_size

	command <<<
		set -o pipefail
		set -e

		cat ${quality_yield} | sed -n '2,2p' > quality_yield.header
		cat ${quality_yield} | sed -n '3,3p' > ${sample}.quality_yield

		cat ${wgs_metrics_algo} | sed -n '2,2p' > wgs_metrics_algo.header
		cat ${wgs_metrics_algo} | sed -n '3,3p' > ${sample}.wgs_metrics_algo

		cat ${aln_metrics} | sed -n '2,2p'  > aln_metrics.header
		cat ${aln_metrics} | sed -n '5,5p'  > ${sample}.aln_metrics

		cat ${is_metrics} | sed -n '2,2p' > is_metrics.header
		cat ${is_metrics} | sed -n '3,3p' > ${sample}.is_metrics

	>>>

	runtime {
		docker:docker
		cluster:cluster_config
		systemDisk:"cloud_ssd 40"
		dataDisk:"cloud_ssd " + disk_size + " /cromwell_root/"
	}

	output {
		File quality_yield_header = "quality_yield.header"
		File quality_yield_data = "${sample}.quality_yield"
		File wgs_metrics_algo_header = "wgs_metrics_algo.header"
		File wgs_metrics_algo_data = "${sample}.wgs_metrics_algo"
		File aln_metrics_header = "aln_metrics.header"
		File aln_metrics_data = "${sample}.aln_metrics"
		File is_metrics_header = "is_metrics.header"
		File is_metrics_data = "${sample}.is_metrics"
	}
}