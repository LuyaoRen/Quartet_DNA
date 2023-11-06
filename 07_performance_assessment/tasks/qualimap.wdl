task qualimap {
	File bam
	File bai
	String bamname = basename(bam,".bam")
	String docker
	String cluster_config
	String disk_size

	command <<<
		set -o pipefail
		set -e
		nt=$(nproc)

		/opt/qualimap/qualimap bamqc -bam ${bam} -outformat PDF:HTML -nt $nt -outdir ${bamname} --java-mem-size=60G

		tar -zcvf ${bamname}_qualimap.zip ${bamname}
	>>>

	runtime {
		docker:docker
		cluster:cluster_config
		systemDisk:"cloud_ssd 40"
		dataDisk:"cloud_ssd " + disk_size + " /cromwell_root/"
	}

	output {
		File zip = "${bamname}_qualimap.zip"
	}
}