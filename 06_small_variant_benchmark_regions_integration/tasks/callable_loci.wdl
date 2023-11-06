task callable_loci {
	
	File LCL5_callable_bed
	File LCL6_callable_bed
	File LCL7_callable_bed
	File LCL8_callable_bed
	String docker
	String disk_size
	String cluster_config
	

	command <<<
		
		cat ${LCL5_callable_bed} | grep -w '^chr1\|^chr2\|^chr3\|^chr4\|^chr5\|^chr6\|^chr7\|^chr8\|^chr9\|^chr10\|^chr11\|^chr12\|^chr13\|^chr14\|^chr15\|^chr16\|^chr17\|^chr18\|^chr19\|^chr20\|^chr21\|^chr22\|^chrX' > LCL5_callable_bed.chr1-22.x.bed


		cat ${LCL6_callable_bed} | grep -w '^chr1\|^chr2\|^chr3\|^chr4\|^chr5\|^chr6\|^chr7\|^chr8\|^chr9\|^chr10\|^chr11\|^chr12\|^chr13\|^chr14\|^chr15\|^chr16\|^chr17\|^chr18\|^chr19\|^chr20\|^chr21\|^chr22\|^chrX' > LCL6_callable_bed.chr1-22.x.bed

		cat ${LCL7_callable_bed} | grep -w '^chr1\|^chr2\|^chr3\|^chr4\|^chr5\|^chr6\|^chr7\|^chr8\|^chr9\|^chr10\|^chr11\|^chr12\|^chr13\|^chr14\|^chr15\|^chr16\|^chr17\|^chr18\|^chr19\|^chr20\|^chr21\|^chr22\|^chrX' > LCL7_callable_bed.chr1-22.x.bed

		cat ${LCL8_callable_bed} | grep -w '^chr1\|^chr2\|^chr3\|^chr4\|^chr5\|^chr6\|^chr7\|^chr8\|^chr9\|^chr10\|^chr11\|^chr12\|^chr13\|^chr14\|^chr15\|^chr16\|^chr17\|^chr18\|^chr19\|^chr20\|^chr21\|^chr22\|^chrX' > LCL8_callable_bed.chr1-22.x.bed


		/opt/ccdg/bedtools-2.27.1/bin/bedtools multiinter -i LCL5_callable_bed.chr1-22.x.bed LCL6_callable_bed.chr1-22.x.bed LCL7_callable_bed.chr1-22.x.bed LCL8_callable_bed.chr1-22.x.bed > Quartet.callable.merged.bed

		cat Quartet.callable.merged.bed | grep "1,2,3,4" | cut -f1-3 > Quartet.callable.merged.intersect.bed

		>>>

		runtime {
		docker:docker
    	cluster:cluster_config
    	systemDisk: "cloud_ssd 40"
    	dataDisk: "cloud_ssd " + disk_size + " /cromwell_root/"
	}

	output {
		File LCL5_callable_filtered_bed = "LCL5_callable_bed.chr1-22.x.bed"
		File LCL6_callable_filtered_bed = "LCL6_callable_bed.chr1-22.x.bed"
		File LCL7_callable_filtered_bed = "LCL7_callable_bed.chr1-22.x.bed"
		File LCL8_callable_filtered_bed = "LCL8_callable_bed.chr1-22.x.bed"
		File callable_merged_bed = "Quartet.callable.merged.bed"
		File callable_merged_intersect_bed = "Quartet.callable.merged.intersect.bed"
	}
}




