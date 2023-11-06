import "./tasks/mapping.wdl" as mapping
import "./tasks/Metrics.wdl" as Metrics
import "./tasks/Dedup.wdl" as Dedup
import "./tasks/deduped_Metrics.wdl" as deduped_Metrics
import "./tasks/Realigner.wdl" as Realigner
import "./tasks/BQSR.wdl" as BQSR
import "./tasks/Haplotyper.wdl" as Haplotyper



workflow {{ project_name }} {

	File inputSamplesFile
	Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)

	String SENTIEON_INSTALL_DIR
	String docker
	
	String fasta
	File ref_dir
	File dbmills_dir
	String db_mills
	File dbsnp_dir
	String dbsnp
	String disk_size
	String cluster_config

	scatter (quartet in inputSamples){

		call mapping.mapping as mapping {
			input: 
			SENTIEON_INSTALL_DIR=SENTIEON_INSTALL_DIR,
			group=quartet[2],
			sample=quartet[2],
			pl="ILLUMINAL",
			fasta=fasta,
			ref_dir=ref_dir,
			fastq_1=quartet[0],
			fastq_2=quartet[1],
			docker=docker,
			disk_size=disk_size,
			cluster_config=cluster_config
		}

		call Metrics.Metrics as Metrics {
			input:
			SENTIEON_INSTALL_DIR=SENTIEON_INSTALL_DIR,
			fasta=fasta,
			ref_dir=ref_dir,
			sorted_bam=mapping.sorted_bam,
			sorted_bam_index=mapping.sorted_bam_index,		
			sample=quartet[2],
			docker=docker,
			disk_size=disk_size,
			cluster_config=cluster_config
		}

		call Dedup.Dedup as Dedup {
			input:
			SENTIEON_INSTALL_DIR=SENTIEON_INSTALL_DIR,
			sorted_bam=mapping.sorted_bam,
			sorted_bam_index=mapping.sorted_bam_index,
			sample=quartet[2],
			docker=docker,
			disk_size=disk_size,
			cluster_config=cluster_config
		}
		call deduped_Metrics.deduped_Metrics as deduped_Metrics {
			input:
			SENTIEON_INSTALL_DIR=SENTIEON_INSTALL_DIR,
			fasta=fasta,
			ref_dir=ref_dir,
			Dedup_bam=Dedup.Dedup_bam,
			Dedup_bam_index=Dedup.Dedup_bam_index,
			sample=quartet[2],
			docker=docker,
			disk_size=disk_size,
			cluster_config=cluster_config
		}
		call Realigner.Realigner as Realigner {
			input:
			SENTIEON_INSTALL_DIR=SENTIEON_INSTALL_DIR,
			fasta=fasta,
			ref_dir=ref_dir,
			Dedup_bam=Dedup.Dedup_bam,
			Dedup_bam_index=Dedup.Dedup_bam_index,
			db_mills=db_mills,
			dbmills_dir=dbmills_dir,
			sample=quartet[2],
			docker=docker,
			disk_size=disk_size,
			cluster_config=cluster_config
		}

		call BQSR.BQSR as BQSR {
			input:
			SENTIEON_INSTALL_DIR=SENTIEON_INSTALL_DIR,
			fasta=fasta,
			ref_dir=ref_dir,
			realigned_bam=Realigner.realigner_bam,
			realigned_bam_index=Realigner.realigner_bam_index,
			db_mills=db_mills,
			dbmills_dir=dbmills_dir,
			dbsnp=dbsnp,
			dbsnp_dir=dbsnp_dir,
			sample=quartet[2],
			docker=docker,
			disk_size=disk_size,
			cluster_config=cluster_config
		}
		call Haplotyper.Haplotyper as Haplotyper {
			input:
			SENTIEON_INSTALL_DIR=SENTIEON_INSTALL_DIR,
			fasta=fasta,
			ref_dir=ref_dir,
			recaled_bam=BQSR.recaled_bam,
			recaled_bam_index=BQSR.recaled_bam_index,
			sample=quartet[2],
			docker=docker,
			disk_size=disk_size,
			cluster_config=cluster_config
		}
	}
}