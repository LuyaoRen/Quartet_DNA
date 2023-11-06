task Metrics {


    File ref_dir
	String SENTIEON_INSTALL_DIR
	String sample
	String docker
	String cluster_config

	String fasta
	File sorted_bam
	File sorted_bam_index
	String disk_size
	


	command <<<
		set -o pipefail
		set -e
		export SENTIEON_LICENSE=192.168.0.55:8990
		nt=$(nproc)
		${SENTIEON_INSTALL_DIR}/bin/sentieon driver -r ${ref_dir}/${fasta} -t $nt -i ${sorted_bam} --algo MeanQualityByCycle ${sample}_mq_metrics.txt --algo QualDistribution ${sample}_qd_metrics.txt --algo GCBias --summary ${sample}_gc_summary.txt ${sample}_gc_metrics.txt --algo AlignmentStat ${sample}_aln_metrics.txt --algo InsertSizeMetricAlgo ${sample}_is_metrics.txt --algo CoverageMetrics --omit_base_output ${sample}_coverage_metrics

		${SENTIEON_INSTALL_DIR}/bin/sentieon plot metrics -o ${sample}_metrics_report.pdf gc=${sample}_gc_metrics.txt qd=${sample}_qd_metrics.txt mq=${sample}_mq_metrics.txt isize=${sample}_is_metrics.txt
	>>>
	
	runtime {
		docker:docker
    	cluster: cluster_config
    	systemDisk: "cloud_ssd 40"
    	dataDisk: "cloud_ssd " + disk_size + " /cromwell_root/"
	}
	output {
		File qd_metrics = "${sample}_qd_metrics.txt"
		File qd_metrics_pdf = "${sample}_qd_metrics.pdf"
		File mq_metrics = "${sample}_mq_metrics.txt"
		File mq_metrics_pdf = "${sample}_mq_metrics.pdf"
		File is_metrics = "${sample}_is_metrics.txt"
		File is_metrics_pdf = "${sample}_is_metrics.pdf"
		File gc_summary = "${sample}_gc_summary.txt"
		File gc_metrics = "${sample}_gc_metrics.txt"
		File gc_metrics_pdf = "${sample}_gc_metrics.pdf"
		File aln_metrics = "${sample}_aln_metrics.txt"
		File coverage_metrics_sample_summary = "${sample}_coverage_metrics.sample_summary"
		File coverage_metrics_sample_statistics = "${sample}_coverage_metrics.sample_statistics"
		File coverage_metrics_sample_interval_statistics = "${sample}_coverage_metrics.sample_interval_statistics"
		File coverage_metrics_sample_cumulative_coverage_proportions = "${sample}_coverage_metrics.sample_cumulative_coverage_proportions"
		File coverage_metrics_sample_cumulative_coverage_counts = "${sample}_coverage_metrics.sample_cumulative_coverage_counts"
	}

}





