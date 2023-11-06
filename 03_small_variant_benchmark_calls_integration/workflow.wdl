import "./tasks/replicate_consensus.wdl" as replicate_consensus
import "./tasks/mendelian.wdl" as mendelian
import "./tasks/filtered_benchmark.wdl" as filtered_benchmark
import "./tasks/merge_mendelian.wdl" as merge_mendelian

workflow {{ project_name }} {
	File ref_dir
	File LCL5_info_dir
	File LCL5_benchmark_call
	File chromo_file
	String fasta
	String chromo
	String cluster_config
	String disk_size


	call replicate_consensus.replicate_consensus as replicate_consensus {
		input:
		chromo_file=chromo_file,
		chromo=chromo,
		cluster_config=cluster_config,
		disk_size=disk_size
	}
	call mendelian.mendelian as mendelian {
		input:
		consensus_vcf=replicate_consensus.consensus_vcf,
		ref_dir=ref_dir,
		fasta=fasta,
		chromo=chromo,
		cluster_config=cluster_config,
		disk_size=disk_size
	}
	call filtered_benchmark.filtered_benchmark as filtered_benchmark {
		input:
		chromo_consensus=replicate_consensus.chromo_consensus,
		LCL5_info_dir=LCL5_info_dir,
		LCL5_benchmark_call=LCL5_benchmark_call,
		chromo=chromo,
		cluster_config=cluster_config,
		disk_size=disk_size
	}
	call merge_mendelian.merge_mendelian as merge_mendelian {
		input:
		D5_trio_vcf=mendelian.D5_trio_vcf,
		D6_trio_vcf=mendelian.D6_trio_vcf,
		consensus_vcf=replicate_consensus.consensus_vcf,
		chromo=chromo,
		cluster_config=cluster_config,
		disk_size=disk_size
	}
}