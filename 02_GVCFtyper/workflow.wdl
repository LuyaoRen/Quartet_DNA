import "./tasks/gVCF_chromo_split.wdl" as gVCF_chromo_split
import "./tasks/GVCFtyper.wdl" as GVCFtyper

workflow {{ project_name }} {

	File inputSamplesFile
	Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)

	String SENTIEON_INSTALL_DIR
	String docker
	String project
	String fasta
	File ref_dir
	String disk_size
	String cluster_config

	scatter (quartet in inputSamples){
		call gVCF_chromo_split.gVCF_chromo_split as gVCF_chromo_split {
			input:
			gvcf=quartet[0],
			gvcf_idx=quartet[1],
			sample=quartet[2],
			SENTIEON_INSTALL_DIR=SENTIEON_INSTALL_DIR,
			docker=docker,
			disk_size=disk_size,
			cluster_config=cluster_config
			}
	}
	call GVCFtyper.GVCFtyper as GVCFtyper_chr1 {
		input:
		ref_dir=ref_dir,
		SENTIEON_INSTALL_DIR=SENTIEON_INSTALL_DIR,
		fasta=fasta,
		vcf=gVCF_chromo_split.chr1_gvcf,
		vcf_idx=gVCF_chromo_split.chr1_gvcf_idx,
		project=project,
		docker=docker,
		cluster_config=cluster_config,
		disk_size=disk_size
	}
	call GVCFtyper.GVCFtyper as GVCFtyper_chr2 {
		input:
		ref_dir=ref_dir,
		SENTIEON_INSTALL_DIR=SENTIEON_INSTALL_DIR,
		fasta=fasta,
		vcf=gVCF_chromo_split.chr2_gvcf,
		vcf_idx=gVCF_chromo_split.chr2_gvcf_idx,
		project=project,
		docker=docker,
		cluster_config=cluster_config,
		disk_size=disk_size
	}
	call GVCFtyper.GVCFtyper as GVCFtyper_chr3 {
		input:
		ref_dir=ref_dir,
		SENTIEON_INSTALL_DIR=SENTIEON_INSTALL_DIR,
		fasta=fasta,
		vcf=gVCF_chromo_split.chr3_gvcf,
		vcf_idx=gVCF_chromo_split.chr3_gvcf_idx,
		project=project,
		docker=docker,
		cluster_config=cluster_config,
		disk_size=disk_size
	}
	call GVCFtyper.GVCFtyper as GVCFtyper_chr4 {
		input:
		ref_dir=ref_dir,
		SENTIEON_INSTALL_DIR=SENTIEON_INSTALL_DIR,
		fasta=fasta,
		vcf=gVCF_chromo_split.chr4_gvcf,
		vcf_idx=gVCF_chromo_split.chr4_gvcf_idx,
		project=project,
		docker=docker,
		cluster_config=cluster_config,
		disk_size=disk_size
	}
	call GVCFtyper.GVCFtyper as GVCFtyper_chr5 {
		input:
		ref_dir=ref_dir,
		SENTIEON_INSTALL_DIR=SENTIEON_INSTALL_DIR,
		fasta=fasta,
		vcf=gVCF_chromo_split.chr5_gvcf,
		vcf_idx=gVCF_chromo_split.chr5_gvcf_idx,
		project=project,
		docker=docker,
		cluster_config=cluster_config,
		disk_size=disk_size
	}
	call GVCFtyper.GVCFtyper as GVCFtyper_chr6 {
		input:
		ref_dir=ref_dir,
		SENTIEON_INSTALL_DIR=SENTIEON_INSTALL_DIR,
		fasta=fasta,
		vcf=gVCF_chromo_split.chr6_gvcf,
		vcf_idx=gVCF_chromo_split.chr6_gvcf_idx,
		project=project,
		docker=docker,
		cluster_config=cluster_config,
		disk_size=disk_size
	}
	call GVCFtyper.GVCFtyper as GVCFtyper_chr7 {
		input:
		ref_dir=ref_dir,
		SENTIEON_INSTALL_DIR=SENTIEON_INSTALL_DIR,
		fasta=fasta,
		vcf=gVCF_chromo_split.chr7_gvcf,
		vcf_idx=gVCF_chromo_split.chr7_gvcf_idx,
		project=project,
		docker=docker,
		cluster_config=cluster_config,
		disk_size=disk_size
	}
	call GVCFtyper.GVCFtyper as GVCFtyper_chr8 {
		input:
		ref_dir=ref_dir,
		SENTIEON_INSTALL_DIR=SENTIEON_INSTALL_DIR,
		fasta=fasta,
		vcf=gVCF_chromo_split.chr8_gvcf,
		vcf_idx=gVCF_chromo_split.chr8_gvcf_idx,
		project=project,
		docker=docker,
		cluster_config=cluster_config,
		disk_size=disk_size
	}
	call GVCFtyper.GVCFtyper as GVCFtyper_chr9 {
		input:
		ref_dir=ref_dir,
		SENTIEON_INSTALL_DIR=SENTIEON_INSTALL_DIR,
		fasta=fasta,
		vcf=gVCF_chromo_split.chr9_gvcf,
		vcf_idx=gVCF_chromo_split.chr9_gvcf_idx,
		project=project,
		docker=docker,
		cluster_config=cluster_config,
		disk_size=disk_size
	}
	call GVCFtyper.GVCFtyper as GVCFtyper_chr10 {
		input:
		ref_dir=ref_dir,
		SENTIEON_INSTALL_DIR=SENTIEON_INSTALL_DIR,
		fasta=fasta,
		vcf=gVCF_chromo_split.chr10_gvcf,
		vcf_idx=gVCF_chromo_split.chr10_gvcf_idx,
		project=project,
		docker=docker,
		cluster_config=cluster_config,
		disk_size=disk_size
	}
	call GVCFtyper.GVCFtyper as GVCFtyper_chr11 {
		input:
		ref_dir=ref_dir,
		SENTIEON_INSTALL_DIR=SENTIEON_INSTALL_DIR,
		fasta=fasta,
		vcf=gVCF_chromo_split.chr11_gvcf,
		vcf_idx=gVCF_chromo_split.chr11_gvcf_idx,
		project=project,
		docker=docker,
		cluster_config=cluster_config,
		disk_size=disk_size
	}
	call GVCFtyper.GVCFtyper as GVCFtyper_chr12 {
		input:
		ref_dir=ref_dir,
		SENTIEON_INSTALL_DIR=SENTIEON_INSTALL_DIR,
		fasta=fasta,
		vcf=gVCF_chromo_split.chr12_gvcf,
		vcf_idx=gVCF_chromo_split.chr12_gvcf_idx,
		project=project,
		docker=docker,
		cluster_config=cluster_config,
		disk_size=disk_size
	}
	call GVCFtyper.GVCFtyper as GVCFtyper_chr13 {
		input:
		ref_dir=ref_dir,
		SENTIEON_INSTALL_DIR=SENTIEON_INSTALL_DIR,
		fasta=fasta,
		vcf=gVCF_chromo_split.chr13_gvcf,
		vcf_idx=gVCF_chromo_split.chr13_gvcf_idx,
		project=project,
		docker=docker,
		cluster_config=cluster_config,
		disk_size=disk_size
	}
	call GVCFtyper.GVCFtyper as GVCFtyper_chr14 {
		input:
		ref_dir=ref_dir,
		SENTIEON_INSTALL_DIR=SENTIEON_INSTALL_DIR,
		fasta=fasta,
		vcf=gVCF_chromo_split.chr14_gvcf,
		vcf_idx=gVCF_chromo_split.chr14_gvcf_idx,
		project=project,
		docker=docker,
		cluster_config=cluster_config,
		disk_size=disk_size
	}
	call GVCFtyper.GVCFtyper as GVCFtyper_chr15 {
		input:
		ref_dir=ref_dir,
		SENTIEON_INSTALL_DIR=SENTIEON_INSTALL_DIR,
		fasta=fasta,
		vcf=gVCF_chromo_split.chr15_gvcf,
		vcf_idx=gVCF_chromo_split.chr15_gvcf_idx,
		project=project,
		docker=docker,
		cluster_config=cluster_config,
		disk_size=disk_size
	}
	call GVCFtyper.GVCFtyper as GVCFtyper_chr16 {
		input:
		ref_dir=ref_dir,
		SENTIEON_INSTALL_DIR=SENTIEON_INSTALL_DIR,
		fasta=fasta,
		vcf=gVCF_chromo_split.chr16_gvcf,
		vcf_idx=gVCF_chromo_split.chr16_gvcf_idx,
		project=project,
		docker=docker,
		cluster_config=cluster_config,
		disk_size=disk_size
	}
	call GVCFtyper.GVCFtyper as GVCFtyper_chr17 {
		input:
		ref_dir=ref_dir,
		SENTIEON_INSTALL_DIR=SENTIEON_INSTALL_DIR,
		fasta=fasta,
		vcf=gVCF_chromo_split.chr17_gvcf,
		vcf_idx=gVCF_chromo_split.chr17_gvcf_idx,
		project=project,
		docker=docker,
		cluster_config=cluster_config,
		disk_size=disk_size
	}
	call GVCFtyper.GVCFtyper as GVCFtyper_chr18 {
		input:
		ref_dir=ref_dir,
		SENTIEON_INSTALL_DIR=SENTIEON_INSTALL_DIR,
		fasta=fasta,
		vcf=gVCF_chromo_split.chr18_gvcf,
		vcf_idx=gVCF_chromo_split.chr18_gvcf_idx,
		project=project,
		docker=docker,
		cluster_config=cluster_config,
		disk_size=disk_size
	}
	call GVCFtyper.GVCFtyper as GVCFtyper_chr19 {
		input:
		ref_dir=ref_dir,
		SENTIEON_INSTALL_DIR=SENTIEON_INSTALL_DIR,
		fasta=fasta,
		vcf=gVCF_chromo_split.chr19_gvcf,
		vcf_idx=gVCF_chromo_split.chr19_gvcf_idx,
		project=project,
		docker=docker,
		cluster_config=cluster_config,
		disk_size=disk_size
	}
	call GVCFtyper.GVCFtyper as GVCFtyper_chr20 {
		input:
		ref_dir=ref_dir,
		SENTIEON_INSTALL_DIR=SENTIEON_INSTALL_DIR,
		fasta=fasta,
		vcf=gVCF_chromo_split.chr20_gvcf,
		vcf_idx=gVCF_chromo_split.chr20_gvcf_idx,
		project=project,
		docker=docker,
		cluster_config=cluster_config,
		disk_size=disk_size
	}
	call GVCFtyper.GVCFtyper as GVCFtyper_chr21 {
		input:
		ref_dir=ref_dir,
		SENTIEON_INSTALL_DIR=SENTIEON_INSTALL_DIR,
		fasta=fasta,
		vcf=gVCF_chromo_split.chr21_gvcf,
		vcf_idx=gVCF_chromo_split.chr21_gvcf_idx,
		project=project,
		docker=docker,
		cluster_config=cluster_config,
		disk_size=disk_size
	}
	call GVCFtyper.GVCFtyper as GVCFtyper_chr22 {
		input:
		ref_dir=ref_dir,
		SENTIEON_INSTALL_DIR=SENTIEON_INSTALL_DIR,
		fasta=fasta,
		vcf=gVCF_chromo_split.chr22_gvcf,
		vcf_idx=gVCF_chromo_split.chr22_gvcf_idx,
		project=project,
		docker=docker,
		cluster_config=cluster_config,
		disk_size=disk_size
	}
	call GVCFtyper.GVCFtyper as GVCFtyper_chrX {
		input:
		ref_dir=ref_dir,
		SENTIEON_INSTALL_DIR=SENTIEON_INSTALL_DIR,
		fasta=fasta,
		vcf=gVCF_chromo_split.chrX_gvcf,
		vcf_idx=gVCF_chromo_split.chrX_gvcf_idx,
		project=project,
		docker=docker,
		cluster_config=cluster_config,
		disk_size=disk_size
	}
}