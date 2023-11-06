task benchmark {
	File filtered_vcf
	File benchmarking_dir
	File ref_dir
	String sample = basename(filtered_vcf,".filtered.vcf")
	String fasta
	String docker
	String cluster_config
	String disk_size


	command <<<
		set -o pipefail
		set -e
		nt=$(nproc)
		mkdir -p /cromwell_root/tmp
		cp -r ${ref_dir} /cromwell_root/tmp/
		cp -r ${benchmarking_dir} /cromwell_root/tmp/

		export HGREF=/cromwell_root/tmp/reference_data/GRCh38.d1.vd1.fa


		echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tLCL5" > LCL5_name
		echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tLCL6" > LCL6_name
		echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tLCL7" > LCL7_name
		echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tLCL8" > LCL8_name


		if [[ ${sample} =~ "LCL5" ]];then
			/opt/hap.py/bin/hap.py /cromwell_root/tmp/reference_datasets_v202103/LCL5.high.confidence.calls.vcf ${filtered_vcf} -f /cromwell_root/tmp/reference_datasets_v202103/Quartet.high.confidence.region.v202103.bed --threads $nt -o ${sample} -r ${ref_dir}/${fasta}
			cat ${filtered_vcf} | grep '##' > header
			cat ${filtered_vcf} | grep -v '#' > body
			cat header LCL5_name body > LCL5.vcf
			/opt/rtg-tools/dist/rtg-tools-3.10.1-4d58ead/rtg bgzip LCL5.vcf -c > ${sample}.reformed.vcf.gz
			/opt/rtg-tools/dist/rtg-tools-3.10.1-4d58ead/rtg index -f vcf ${sample}.reformed.vcf.gz
		elif [[ ${sample} =~ "LCL6" ]]; then
		    /opt/hap.py/bin/hap.py /cromwell_root/tmp/reference_datasets_v202103/LCL6.high.confidence.calls.vcf ${filtered_vcf} -f /cromwell_root/tmp/reference_datasets_v202103/Quartet.high.confidence.region.v202103.bed --threads $nt -o ${sample} -r ${ref_dir}/${fasta}
			cat ${filtered_vcf} | grep '##' > header
			cat ${filtered_vcf} | grep -v '#' > body
			cat header LCL6_name body > LCL6.vcf
			/opt/rtg-tools/dist/rtg-tools-3.10.1-4d58ead/rtg bgzip LCL6.vcf -c > ${sample}.reformed.vcf.gz
			/opt/rtg-tools/dist/rtg-tools-3.10.1-4d58ead/rtg index -f vcf ${sample}.reformed.vcf.gz
	    elif [[ ${sample} =~ "LCL7" ]]; then
	        /opt/hap.py/bin/hap.py /cromwell_root/tmp/reference_datasets_v202103/LCL7.high.confidence.calls.vcf ${filtered_vcf} -f /cromwell_root/tmp/reference_datasets_v202103/Quartet.high.confidence.region.v202103.bed --threads $nt -o ${sample} -r ${ref_dir}/${fasta}
			cat ${filtered_vcf} | grep '##' > header
			cat ${filtered_vcf} | grep -v '#' > body
			cat header LCL7_name body > LCL7.vcf
			/opt/rtg-tools/dist/rtg-tools-3.10.1-4d58ead/rtg bgzip LCL7.vcf -c > ${sample}.reformed.vcf.gz
			/opt/rtg-tools/dist/rtg-tools-3.10.1-4d58ead/rtg index -f vcf ${sample}.reformed.vcf.gz
		elif [[ ${sample} =~ "LCL8" ]]; then
			/opt/hap.py/bin/hap.py /cromwell_root/tmp/reference_datasets_v202103/LCL8.high.confidence.calls.vcf ${filtered_vcf} -f /cromwell_root/tmp/reference_datasets_v202103/Quartet.high.confidence.region.v202103.bed --threads $nt -o ${sample} -r ${ref_dir}/${fasta}
			cat ${filtered_vcf} | grep '##' > header
			cat ${filtered_vcf} | grep -v '#' > body
			cat header LCL8_name body > LCL8.vcf
			/opt/rtg-tools/dist/rtg-tools-3.10.1-4d58ead/rtg bgzip LCL8.vcf -c > ${sample}.reformed.vcf.gz
			/opt/rtg-tools/dist/rtg-tools-3.10.1-4d58ead/rtg index -f vcf ${sample}.reformed.vcf.gz
	    else
	        echo "only for quartet samples"
	    fi
	>>>

	runtime {
		docker:docker
		cluster:cluster_config
		systemDisk:"cloud_ssd 40"
		dataDisk:"cloud_ssd " + disk_size + " /cromwell_root/"
	}

	output {
		File rtg_vcf = "${sample}.reformed.vcf.gz"
		File rtg_vcf_index = "${sample}.reformed.vcf.gz.tbi"
		File gzip_vcf = "${sample}.vcf.gz"
		File gzip_vcf_index = "${sample}.vcf.gz.tbi"
		File roc_all_csv = "${sample}.roc.all.csv.gz"
		File roc_indel = "${sample}.roc.Locations.INDEL.csv.gz"
		File roc_indel_pass = "${sample}.roc.Locations.INDEL.PASS.csv.gz"
		File roc_snp = "${sample}.roc.Locations.SNP.csv.gz"
		File roc_snp_pass = "${sample}.roc.Locations.SNP.PASS.csv.gz"
		File summary = "${sample}.summary.csv"
		File extended = "${sample}.extended.csv"
		File metrics = "${sample}.metrics.json.gz"
	}
}