task merge_mendelian {
	File D5_trio_vcf
	File D6_trio_vcf
	File family_vcf
	String family_name = basename(family_vcf,".family.vcf")
	String docker
	String cluster_config
	String disk_size
	
	command <<<
		cat ${D5_trio_vcf} | grep -v '##' > ${family_name}.D5.txt
		cat ${D6_trio_vcf} | grep -v '##' > ${family_name}.D6.txt
		cat ${family_vcf} | grep -v '##' | awk '
		    BEGIN { OFS = "\t" }
		    NF > 2 && FNR > 1 { 
		        for ( i=9; i<=NF; i++ ) { 
		            split($i,a,":") ;$i = a[1];
		        } 
		    } 
		    { print }
		' > ${family_name}.consensus.txt
		python /opt/merge_two_family_with_genotype.py -LCL5 ${family_name}.D5.txt -LCL6 ${family_name}.D6.txt -genotype ${family_name}.consensus.txt -family ${family_name}
	>>>

	runtime {
		docker:docker
		cluster: cluster_config
		systemDisk: "cloud_ssd 40"
		dataDisk: "cloud_ssd " + disk_size + " /cromwell_root/"
	}
	output {
		File project_mendelian = "${family_name}.txt"
		File project_mendelian_summary = "${family_name}.summary.txt"
	}
}