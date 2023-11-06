task BQSR {
	
    File ref_dir
    File dbsnp_dir
    File dbmills_dir
	String SENTIEON_INSTALL_DIR
	String fasta
	String dbsnp
	String db_mills
	File realigned_bam
	File realigned_bam_index
	String sample = basename(realigned_bam,".sorted.deduped.realigned.bam")
	String docker
	String cluster_config
	String disk_size

	
	command <<<
		set -o pipefail
		set -e
		export SENTIEON_LICENSE=192.168.0.55:8990
		nt=$(nproc)

		${SENTIEON_INSTALL_DIR}/bin/sentieon driver -r ${ref_dir}/${fasta} -t $nt -i ${realigned_bam} --algo QualCal -k ${dbsnp_dir}/${dbsnp} -k ${dbmills_dir}/${db_mills} ${sample}_recal_data.table

		${SENTIEON_INSTALL_DIR}/bin/sentieon driver -r ${ref_dir}/${fasta} -t $nt -i ${realigned_bam} -q ${sample}_recal_data.table --algo QualCal -k ${dbsnp_dir}/${dbsnp} -k ${dbmills_dir}/${db_mills} ${sample}_recal_data.table.post --algo ReadWriter ${sample}.sorted.deduped.realigned.recaled.bam

		${SENTIEON_INSTALL_DIR}/bin/sentieon driver -t $nt --algo QualCal --plot --before ${sample}_recal_data.table --after ${sample}_recal_data.table.post ${sample}_recal_data.csv

		${SENTIEON_INSTALL_DIR}/bin/sentieon plot QualCal -o ${sample}_bqsrreport.pdf ${sample}_recal_data.csv

	>>>
	runtime {
		docker:docker
		cluster: cluster_config
		systemDisk: "cloud_ssd 40"
		dataDisk: "cloud_ssd " + disk_size + " /cromwell_root/"	
	}

	output {
		File recal_table = "${sample}_recal_data.table"
		File recal_post = "${sample}_recal_data.table.post"
		File recaled_bam = "${sample}.sorted.deduped.realigned.recaled.bam"
		File recaled_bam_index = "${sample}.sorted.deduped.realigned.recaled.bam.bai"
		File recal_csv = "${sample}_recal_data.csv"
		File bqsrreport_pdf = "${sample}_bqsrreport.pdf"
	}
}
