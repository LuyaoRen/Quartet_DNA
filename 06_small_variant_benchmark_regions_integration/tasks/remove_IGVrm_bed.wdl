task remove_IGVrm_bed {
	
	File benchmark_region
	File bed_10X
	File PMRA_bed
	String docker
	String disk_size
	String cluster_config
	

	command <<<

		cat ${bed_10X} | cut -f1,13,14 > false.10X.bed

		cat ${PMRA_bed} | cut -f1,14,15 > false.PMRA.bed

		cat false.10X.bed false.PMRA.bed | sort -k1,1 -k2,2n > false.positive.bed

		/opt/ccdg/bedtools-2.27.1/bin/bedtools subtract -a ${benchmark_region} -b false.positive.bed > benchmark_regions.filtered.bed

		>>>

		runtime {
		docker:docker
    	cluster:cluster_config
    	systemDisk: "cloud_ssd 40"
    	dataDisk: "cloud_ssd " + disk_size + " /cromwell_root/"
	}

	output {
		File benchmark_filtered_region = "benchmark_regions.filtered.bed"
	}
}







