task filter_vcf {
	
	File benchmark_region
	File LCL5_vcf
	File LCL6_vcf
	File LCL7_vcf
	File LCL8_vcf
	File LCL5_vcf_idx
	File LCL6_vcf_idx
	File LCL7_vcf_idx
	File LCL8_vcf_idx
	File vcf_info
	String docker
	String disk_size
	String cluster_config
	

	command <<<

		/opt/rtg-tools/dist/rtg-tools-3.10.1-4d58ead/rtg vcffilter -i ${LCL5_vcf} --include-bed=${benchmark_region} -o LCL5.high.confidence.calls.vcf.gz

		/opt/rtg-tools/dist/rtg-tools-3.10.1-4d58ead/rtg vcffilter -i ${LCL6_vcf} --include-bed=${benchmark_region} -o LCL6.high.confidence.calls.vcf.gz

		/opt/rtg-tools/dist/rtg-tools-3.10.1-4d58ead/rtg vcffilter -i ${LCL7_vcf} --include-bed=${benchmark_region} -o LCL7.high.confidence.calls.vcf.gz

		/opt/rtg-tools/dist/rtg-tools-3.10.1-4d58ead/rtg vcffilter -i ${LCL8_vcf} --include-bed=${benchmark_region} -o LCL8.high.confidence.calls.vcf.gz

		zcat LCL5.high.confidence.calls.vcf.gz | grep '#' > LCL5.header
		zcat LCL6.high.confidence.calls.vcf.gz | grep '#' > LCL6.header
		zcat LCL7.high.confidence.calls.vcf.gz | grep '#' > LCL7.header
		zcat LCL8.high.confidence.calls.vcf.gz | grep '#' > LCL8.header

		gunzip LCL5.high.confidence.calls.vcf.gz
		gunzip LCL6.high.confidence.calls.vcf.gz
		gunzip LCL7.high.confidence.calls.vcf.gz
		gunzip LCL8.high.confidence.calls.vcf.gz

		cat LCL5.high.confidence.calls.vcf | grep -v '#' > LCL5.high.confidence.calls.body
		cat LCL6.high.confidence.calls.vcf | grep -v '#' > LCL6.high.confidence.calls.body
		cat LCL7.high.confidence.calls.vcf | grep -v '#' > LCL7.high.confidence.calls.body
		cat LCL8.high.confidence.calls.vcf | grep -v '#' > LCL8.high.confidence.calls.body


		python /opt/annotate_vcf.py -info ${vcf_info} -vcf LCL5.high.confidence.calls.body -prefix LCL5
		python /opt/annotate_vcf.py -info ${vcf_info} -vcf LCL6.high.confidence.calls.body -prefix LCL6
		python /opt/annotate_vcf.py -info ${vcf_info} -vcf LCL7.high.confidence.calls.body -prefix LCL7
		python /opt/annotate_vcf.py -info ${vcf_info} -vcf LCL8.high.confidence.calls.body -prefix LCL8


		cat LCL5.annotated.txt | awk '{print $1"\t"$2"\t.\t"$4"\t"$5"\t.\t.\tVOTE="$13"\tGT:ALT:DP\t"$10":"$18":"$17}' | cat LCL5.header - > LCL5.high.confidence.calls.annotated.vcf
		cat LCL6.annotated.txt | awk '{print $1"\t"$2"\t.\t"$4"\t"$5"\t.\t.\tVOTE="$14"\tGT:ALT:DP\t"$10":"$20":"$19}' | cat LCL6.header - > LCL6.high.confidence.calls.annotated.vcf
		cat LCL7.annotated.txt | awk '{print $1"\t"$2"\t.\t"$4"\t"$5"\t.\t.\tVOTE="$15"\tGT:ALT:DP\t"$10":"$22":"$21}' | cat LCL7.header - > LCL7.high.confidence.calls.annotated.vcf
		cat LCL8.annotated.txt | awk '{print $1"\t"$2"\t.\t"$4"\t"$5"\t.\t.\tVOTE="$16"\tGT:ALT:DP\t"$10":"$24":"$23}' | cat LCL8.header - > LCL8.high.confidence.calls.annotated.vcf


		>>>

		runtime {
		docker:docker
    	cluster:cluster_config
    	systemDisk: "cloud_ssd 40"
    	dataDisk: "cloud_ssd " + disk_size + " /cromwell_root/"
	}

	output {
		File LCL5_filtered_vcf = "LCL5.high.confidence.calls.vcf"
		File LCL6_filtered_vcf = "LCL6.high.confidence.calls.vcf"
		File LCL7_filtered_vcf = "LCL7.high.confidence.calls.vcf"
		File LCL8_filtered_vcf = "LCL8.high.confidence.calls.vcf"
		File LCL5_annotated_vcf = "LCL5.high.confidence.calls.annotated.vcf"
		File LCL6_annotated_vcf = "LCL6.high.confidence.calls.annotated.vcf"
		File LCL7_annotated_vcf = "LCL7.high.confidence.calls.annotated.vcf"
		File LCL8_annotated_vcf = "LCL8.high.confidence.calls.annotated.vcf"

	}
}







