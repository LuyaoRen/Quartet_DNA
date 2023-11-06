task quartet_mendelian {
	Array[File] project_mendelian_summary
	String project
	String docker
	String cluster_config
	String disk_size

	command <<<
		for i in ${sep=" " project_mendelian_summary}
		do
		  cat $i |  sed -n '2,3p' >> mendelian.summary
		done
		sed '1iFamily\tTotal_Variants\tMendelian_Concordant_Variants\tMendelian_Concordance_Rate' mendelian.summary > mendelian.txt

		cat mendelian.txt | grep 'INDEL' | cut -f4 | grep -v 'Mendelian_Concordance_Rate' | awk '{for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} 
          END {for (i=1;i<=NF;i++) {
          printf "%f %f \n", sum[i]/NR, sqrt((sumsq[i]-sum[i]^2/NR)/NR)}
         }' >> quartet_indel_aver-std.txt

		cat mendelian.txt | grep 'SNV' | cut -f4 | grep -v 'Mendelian_Concordance_Rate' | awk '{for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} 
          END {for (i=1;i<=NF;i++) {
          printf "%f %f \n", sum[i]/NR, sqrt((sumsq[i]-sum[i]^2/NR)/NR)}
         }' >> quartet_snv_aver-std.txt


	>>>

	runtime {
		docker:docker
    	cluster:cluster_config
    	systemDisk:"cloud_ssd 40"
    	dataDisk:"cloud_ssd " + disk_size + " /cromwell_root/"
	}
	output {
		File mendelian_summary="mendelian.txt"
		File snv_aver_std = "quartet_snv_aver-std.txt"
		File indel_aver_std = "quartet_indel_aver-std.txt"
	}
}