task remove_IGVrm_vcf {
	
	File benchmark_filtered_region
	File LCL5_annotated_vcf
	File LCL6_annotated_vcf
	File LCL7_annotated_vcf
	File LCL8_annotated_vcf
	String docker
	String disk_size
	String cluster_config
	

	command <<<

		/opt/rtg-tools/dist/rtg-tools-3.10.1-4d58ead/rtg bgzip ${LCL5_annotated_vcf} -c > LCL5.vcf.gz
		/opt/rtg-tools/dist/rtg-tools-3.10.1-4d58ead/rtg index -f vcf LCL5.vcf.gz

		/opt/rtg-tools/dist/rtg-tools-3.10.1-4d58ead/rtg bgzip ${LCL6_annotated_vcf} -c > LCL6.vcf.gz
		/opt/rtg-tools/dist/rtg-tools-3.10.1-4d58ead/rtg index -f vcf LCL6.vcf.gz

		/opt/rtg-tools/dist/rtg-tools-3.10.1-4d58ead/rtg bgzip ${LCL7_annotated_vcf} -c > LCL7.vcf.gz
		/opt/rtg-tools/dist/rtg-tools-3.10.1-4d58ead/rtg index -f vcf LCL7.vcf.gz

		/opt/rtg-tools/dist/rtg-tools-3.10.1-4d58ead/rtg bgzip ${LCL8_annotated_vcf} -c > LCL8.vcf.gz
		/opt/rtg-tools/dist/rtg-tools-3.10.1-4d58ead/rtg index -f vcf LCL8.vcf.gz



		/opt/rtg-tools/dist/rtg-tools-3.10.1-4d58ead/rtg vcffilter -i LCL5.vcf.gz --include-bed=${benchmark_filtered_region} -o LCL5.high.confidence.calls.vcf.gz

		/opt/rtg-tools/dist/rtg-tools-3.10.1-4d58ead/rtg vcffilter -i LCL6.vcf.gz --include-bed=${benchmark_filtered_region} -o LCL6.high.confidence.calls.vcf.gz

		/opt/rtg-tools/dist/rtg-tools-3.10.1-4d58ead/rtg vcffilter -i LCL7.vcf.gz --include-bed=${benchmark_filtered_region} -o LCL7.high.confidence.calls.vcf.gz

		/opt/rtg-tools/dist/rtg-tools-3.10.1-4d58ead/rtg vcffilter -i LCL8.vcf.gz --include-bed=${benchmark_filtered_region} -o LCL8.high.confidence.calls.vcf.gz

		>>>

		runtime {
		docker:docker
    	cluster:cluster_config
    	systemDisk: "cloud_ssd 40"
    	dataDisk: "cloud_ssd " + disk_size + " /cromwell_root/"
	}

	output {
		File LCL5_IGVrm_vcf = "LCL5.high.confidence.calls.vcf.gz"
		File LCL6_IGVrm_vcf = "LCL6.high.confidence.calls.vcf.gz"
		File LCL7_IGVrm_vcf = "LCL7.high.confidence.calls.vcf.gz"
		File LCL8_IGVrm_vcf = "LCL8.high.confidence.calls.vcf.gz"
	}
}







