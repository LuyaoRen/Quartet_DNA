task final_merge {
	
	File callable_merged_intersect_bed
	File HR_merged_intersect_bed
	File variants_merged_bed
	String docker
	String disk_size
	String cluster_config
	

	command <<<

		cat ${HR_merged_intersect_bed} ${variants_merged_bed} | sort -k1,1 -k2,2n | /opt/ccdg/bedtools-2.27.1/bin/bedtools merge -i - > variant_invariant.bed

		/opt/ccdg/bedtools-2.27.1/bin/bedtools intersect -a variant_invariant.bed -b ${callable_merged_intersect_bed} > benchmark_regions.bed

		>>>

		runtime {
		docker:docker
    	cluster:cluster_config
    	systemDisk: "cloud_ssd 40"
    	dataDisk: "cloud_ssd " + disk_size + " /cromwell_root/"
	}

	output {
		File variant_invariant = "variant_invariant.bed"
		File benchmark_region = "benchmark_regions.bed"
	}
}







