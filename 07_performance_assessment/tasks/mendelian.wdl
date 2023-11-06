task mendelian {
	File family_vcf
	File ref_dir
	String family_name = basename(family_vcf,".family.vcf")
	String fasta
	String docker
	String cluster_config
	String disk_size
	
	command <<<
		export LD_LIBRARY_PATH=/opt/htslib-1.9
		nt=$(nproc)

		echo -e "${family_name}\tLCL8\t0\t0\t2\t-9\n${family_name}\tLCL7\t0\t0\t1\t-9\n${family_name}\tLCL5\tLCL7\tLCL8\t2\t-9" > ${family_name}.D5.ped

		mkdir VBT_D5
		/opt/VBT-TrioAnalysis/vbt mendelian -ref ${ref_dir}/${fasta} -mother ${family_vcf} -father ${family_vcf} -child ${family_vcf} -pedigree ${family_name}.D5.ped -outDir VBT_D5 -out-prefix ${family_name}.D5 --output-violation-regions -thread-count $nt

		cat VBT_D5/${family_name}.D5_trio.vcf > ${family_name}.D5.vcf

		echo -e "${family_name}\tLCL8\t0\t0\t2\t-9\n${family_name}\tLCL7\t0\t0\t1\t-9\n${family_name}\tLCL6\tLCL7\tLCL8\t2\t-9" > ${family_name}.D6.ped

		mkdir VBT_D6
		/opt/VBT-TrioAnalysis/vbt mendelian -ref ${ref_dir}/${fasta} -mother ${family_vcf} -father ${family_vcf} -child ${family_vcf} -pedigree ${family_name}.D6.ped -outDir VBT_D6 -out-prefix ${family_name}.D6 --output-violation-regions -thread-count $nt

		cat VBT_D6/${family_name}.D6_trio.vcf > ${family_name}.D6.vcf
	>>>

	runtime {
		docker:docker
		cluster: cluster_config
		systemDisk: "cloud_ssd 40"
		dataDisk: "cloud_ssd " + disk_size + " /cromwell_root/"
	}
	output {
		File D5_ped = "${family_name}.D5.ped"
		File D6_ped = "${family_name}.D6.ped"
		Array[File] D5_mendelian = glob("VBT_D5/*")
		Array[File] D6_mendelian = glob("VBT_D6/*")
		File D5_trio_vcf = "${family_name}.D5.vcf"
		File D6_trio_vcf = "${family_name}.D6.vcf"
	}
}



