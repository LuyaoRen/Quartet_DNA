import "./tasks/CallableLoci.wdl" as CallableLoci
import "./tasks/mergeBed.wdl" as mergeBed
import "./tasks/bedVote.wdl" as bedVote


workflow {{ project_name }} {

	File inputSamplesFile
	Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)
	String quartet_sample
	String disk_size
	String cluster_config
	
	scatter (quartet in inputSamples){
		call CallableLoci.CallableLoci as CallableLoci {
			input: 
			bed=quartet[0],
			sample=quartet[1],
			disk_size=disk_size,
			cluster_config=cluster_config
		}
	}
	call mergeBed.mergeBed as mergeBed {
		input:
		callable_bed=CallableLoci.callable_bed,
		sample=quartet_sample,
		disk_size=disk_size,
		cluster_config=cluster_config
	}
	call bedVote.bedVote as bedVote {
		input:
		merged_bed=mergeBed.merged_bed,
		sample=quartet_sample,
		disk_size=disk_size,
		cluster_config=cluster_config
	}
}

