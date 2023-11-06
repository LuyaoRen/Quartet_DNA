task variant_bed {
	
	File LCL5_variants_bed
	File LCL6_variants_bed
	File LCL7_variants_bed
	File LCL8_variants_bed
	String docker
	String disk_size
	String cluster_config
	

	command <<<
		cat ${LCL5_variants_bed} | cut -f1,11,12 > LCL5_variants_bed.chr1-22.x.bed
		cat ${LCL6_variants_bed} | cut -f1,11,12 > LCL6_variants_bed.chr1-22.x.bed
		cat ${LCL7_variants_bed} | cut -f1,11,12 > LCL7_variants_bed.chr1-22.x.bed
		cat ${LCL8_variants_bed} | cut -f1,11,12 > LCL8_variants_bed.chr1-22.x.bed

		cat LCL5_variants_bed.chr1-22.x.bed LCL6_variants_bed.chr1-22.x.bed LCL7_variants_bed.chr1-22.x.bed LCL8_variants_bed.chr1-22.x.bed | sort -k1,1 -k2,2n | /opt/ccdg/bedtools-2.27.1/bin/bedtools merge -i - > Quartet.variants.merged.union.bed

		>>>

		runtime {
		docker:docker
    	cluster:cluster_config
    	systemDisk: "cloud_ssd 40"
    	dataDisk: "cloud_ssd " + disk_size + " /cromwell_root/"
	}

	output {
		File variants_merged_bed = "Quartet.variants.merged.union.bed"
	}
}







