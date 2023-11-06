task filtered_benchmark {
	File chromo_consensus
	File LCL5_info_dir
	File LCL5_benchmark_call
	String chromo
	String docker
	String cluster_config
	String disk_size
	
	command <<<
		cat ${chromo_consensus} | cut -f1-7 | grep notConsensus > ${chromo}.filtered.txt
		cat ${chromo_consensus} | cut -f1-7 | grep '\./\.' >> ${chromo}.filtered.txt
		zcat ${LCL5_benchmark_call} | grep -v '##' | grep -w ${chromo} > LCL5.benchmark.${chromo}.txt
		
		echo '''Quartet_DNA_BGI_SEQ2000_BGI_LCL5_1_20180518_mendelian_vcfInfo.vcf       Quartet_DNA_BGI_SEQ2000_BGI_LCL5_1_20180518
Quartet_DNA_BGI_SEQ2000_BGI_LCL5_2_20180530_mendelian_vcfInfo.vcf       Quartet_DNA_BGI_SEQ2000_BGI_LCL5_2_20180530
Quartet_DNA_BGI_SEQ2000_BGI_LCL5_3_20180530_mendelian_vcfInfo.vcf       Quartet_DNA_BGI_SEQ2000_BGI_LCL5_3_20180530
Quartet_DNA_BGI_T7_WGE_LCL5_1_20191105_mendelian_vcfInfo.vcf    Quartet_DNA_BGI_T7_WGE_LCL5_1_20191105
Quartet_DNA_BGI_T7_WGE_LCL5_2_20191105_mendelian_vcfInfo.vcf    Quartet_DNA_BGI_T7_WGE_LCL5_2_20191105
Quartet_DNA_BGI_T7_WGE_LCL5_3_20191105_mendelian_vcfInfo.vcf    Quartet_DNA_BGI_T7_WGE_LCL5_3_20191105
Quartet_DNA_ILM_Nova_ARD_LCL5_1_20181108_mendelian_vcfInfo.vcf  Quartet_DNA_ILM_Nova_ARD_LCL5_1_20181108
Quartet_DNA_ILM_Nova_ARD_LCL5_2_20181108_mendelian_vcfInfo.vcf  Quartet_DNA_ILM_Nova_ARD_LCL5_2_20181108
Quartet_DNA_ILM_Nova_ARD_LCL5_3_20181108_mendelian_vcfInfo.vcf  Quartet_DNA_ILM_Nova_ARD_LCL5_3_20181108
Quartet_DNA_ILM_Nova_ARD_LCL5_4_20190111_mendelian_vcfInfo.vcf  Quartet_DNA_ILM_Nova_ARD_LCL5_4_20190111
Quartet_DNA_ILM_Nova_ARD_LCL5_5_20190111_mendelian_vcfInfo.vcf  Quartet_DNA_ILM_Nova_ARD_LCL5_5_20190111
Quartet_DNA_ILM_Nova_ARD_LCL5_6_20190111_mendelian_vcfInfo.vcf  Quartet_DNA_ILM_Nova_ARD_LCL5_6_20190111
Quartet_DNA_ILM_Nova_BRG_LCL5_1_20180930_mendelian_vcfInfo.vcf  Quartet_DNA_ILM_Nova_BRG_LCL5_1_20180930
Quartet_DNA_ILM_Nova_BRG_LCL5_2_20180930_mendelian_vcfInfo.vcf  Quartet_DNA_ILM_Nova_BRG_LCL5_2_20180930
Quartet_DNA_ILM_Nova_BRG_LCL5_3_20180930_mendelian_vcfInfo.vcf  Quartet_DNA_ILM_Nova_BRG_LCL5_3_20180930
Quartet_DNA_ILM_Nova_WUX_LCL5_1_20190917_mendelian_vcfInfo.vcf  Quartet_DNA_ILM_Nova_WUX_LCL5_1_20190917
Quartet_DNA_ILM_Nova_WUX_LCL5_2_20190917_mendelian_vcfInfo.vcf  Quartet_DNA_ILM_Nova_WUX_LCL5_2_20190917
Quartet_DNA_ILM_Nova_WUX_LCL5_3_20190917_mendelian_vcfInfo.vcf  Quartet_DNA_ILM_Nova_WUX_LCL5_3_20190917
Quartet_DNA_ILM_XTen_ARD_LCL5_1_20170403_mendelian_vcfInfo.vcf  Quartet_DNA_ILM_XTen_ARD_LCL5_1_20170403
Quartet_DNA_ILM_XTen_ARD_LCL5_2_20170403_mendelian_vcfInfo.vcf  Quartet_DNA_ILM_XTen_ARD_LCL5_2_20170403
Quartet_DNA_ILM_XTen_ARD_LCL5_3_20170403_mendelian_vcfInfo.vcf  Quartet_DNA_ILM_XTen_ARD_LCL5_3_20170403
Quartet_DNA_ILM_XTen_NVG_LCL5_1_20170329_mendelian_vcfInfo.vcf  Quartet_DNA_ILM_XTen_NVG_LCL5_1_20170329
Quartet_DNA_ILM_XTen_NVG_LCL5_2_20170329_mendelian_vcfInfo.vcf  Quartet_DNA_ILM_XTen_NVG_LCL5_2_20170329
Quartet_DNA_ILM_XTen_NVG_LCL5_3_20170329_mendelian_vcfInfo.vcf  Quartet_DNA_ILM_XTen_NVG_LCL5_3_20170329
Quartet_DNA_ILM_XTen_WUX_LCL5_1_20170216_mendelian_vcfInfo.vcf  Quartet_DNA_ILM_XTen_WUX_LCL5_1_20170216
Quartet_DNA_ILM_XTen_WUX_LCL5_2_20170216_mendelian_vcfInfo.vcf  Quartet_DNA_ILM_XTen_WUX_LCL5_2_20170216
Quartet_DNA_ILM_XTen_WUX_LCL5_3_20170216_mendelian_vcfInfo.vcf  Quartet_DNA_ILM_XTen_WUX_LCL5_3_20170216''' > LCL5_vcf_files

		cat LCL5_vcf_files | while read a b
		do
			cat ${LCL5_info_dir}/$a | grep -v '##' > $b.txt
			python /opt/get_filtered_benchmark_vcfinfo.py -filtered ${chromo}.filtered.txt -benchmark LCL5.benchmark.${chromo}.txt -vcf $b.txt -filename $b.${chromo}
		done

	>>>

	runtime {
		docker:docker
		cluster: cluster_config
		systemDisk: "cloud_ssd 40"
		dataDisk: "cloud_ssd " + disk_size + " /cromwell_root/"
	}

	output {
		File filtered = "${chromo}.filtered.txt"
		File benchmark = "LCL5.benchmark.${chromo}.txt"
		Array[File] filtered_vcf = glob("*.${chromo}.filtered.txt")
		Array[File] benchmark_vcf = glob("*.${chromo}.benchmark.txt")
	}

}