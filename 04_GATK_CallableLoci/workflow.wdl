import "./tasks/CallableLoci.wdl" as CallableLoci

workflow {{ project_name }} {

	File bam
	File bam_index
	File ref_dir
	String fasta
	String sample
	String docker
	String disk_size
	String cluster_config
	
	call CallableLoci.CallableLoci as CallableLoci {
		input: 
		bam=bam,
		bam_index=bam_index,
		ref_dir=ref_dir,
		fasta=fasta,
		sample=sample,
		docker=docker,
		disk_size=disk_size,
		cluster_config=cluster_config
	}

	}

