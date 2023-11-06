import "./tasks/callable_loci.wdl" as callable_loci
import "./tasks/final_merge.wdl" as final_merge
import "./tasks/homo_bed.wdl" as homo_bed
import "./tasks/variant_bed.wdl" as variant_bed
import "./tasks/filter_vcf.wdl" as filter_vcf
import "./tasks/remove_IGVrm_bed.wdl" as remove_IGVrm_bed
import "./tasks/remove_IGVrm_vcf.wdl" as remove_IGVrm_vcf


workflow {{ project_name }} {

	File LCL5_callable_bed
	File LCL6_callable_bed
	File LCL7_callable_bed
	File LCL8_callable_bed
	File LCL5_HR_bed
	File LCL6_HR_bed
	File LCL7_HR_bed
	File LCL8_HR_bed
	File LCL5_variants_bed
	File LCL6_variants_bed
	File LCL7_variants_bed
	File LCL8_variants_bed
	File bed_10X
	File PMRA_bed
	File LCL5_vcf
	File LCL6_vcf
	File LCL7_vcf
	File LCL8_vcf
	File LCL5_vcf_idx
	File LCL6_vcf_idx
	File LCL7_vcf_idx
	File LCL8_vcf_idx
	File vcf_info
	String RTGdocker
	String BEDdocker	
	String disk_size
	String cluster_config
	
	call callable_loci.callable_loci as callable_loci {
		input:
		LCL5_callable_bed=LCL5_callable_bed,
		LCL6_callable_bed=LCL6_callable_bed,
		LCL7_callable_bed=LCL7_callable_bed,
		LCL8_callable_bed=LCL8_callable_bed,
		docker=BEDdocker,
		disk_size=disk_size,
		cluster_config=cluster_config
	}

	call homo_bed.homo_bed as homo_bed {
		input:
		LCL5_HR_bed=LCL5_HR_bed,
		LCL6_HR_bed=LCL6_HR_bed,
		LCL7_HR_bed=LCL7_HR_bed,
		LCL8_HR_bed=LCL8_HR_bed,
		docker=BEDdocker,
		disk_size=disk_size,
		cluster_config=cluster_config
	}

	call variant_bed.variant_bed as variant_bed {
		input:
		LCL5_variants_bed=LCL5_variants_bed,
		LCL6_variants_bed=LCL6_variants_bed,
		LCL7_variants_bed=LCL7_variants_bed,
		LCL8_variants_bed=LCL8_variants_bed,
		docker=BEDdocker,
		disk_size=disk_size,
		cluster_config=cluster_config
	}

	call final_merge.final_merge as final_merge {
		input:
		callable_merged_intersect_bed=callable_loci.callable_merged_intersect_bed,
		HR_merged_intersect_bed=homo_bed.HR_merged_intersect_bed,
		variants_merged_bed=variant_bed.variants_merged_bed,
		docker=BEDdocker,
		disk_size=disk_size,
		cluster_config=cluster_config
	}

	call filter_vcf.filter_vcf as filter_vcf {
		input:
		benchmark_region=final_merge.benchmark_region,
		vcf_info=vcf_info,
		LCL5_vcf=LCL5_vcf,
		LCL6_vcf=LCL6_vcf,
		LCL7_vcf=LCL7_vcf,
		LCL8_vcf=LCL8_vcf,
		LCL5_vcf_idx=LCL5_vcf_idx,
		LCL6_vcf_idx=LCL6_vcf_idx,	
		LCL7_vcf_idx=LCL7_vcf_idx,
		LCL8_vcf_idx=LCL8_vcf_idx,
		docker=RTGdocker,
		disk_size=disk_size,
		cluster_config=cluster_config
	}

	call remove_IGVrm_bed.remove_IGVrm_bed as remove_IGVrm_bed{
		input:
		benchmark_region=final_merge.benchmark_region,
		bed_10X=bed_10X,
		PMRA_bed=PMRA_bed,
		docker=BEDdocker,
		disk_size=disk_size,
		cluster_config=cluster_config
	}

	call remove_IGVrm_vcf.remove_IGVrm_vcf as remove_IGVrm_vcf{
		input:	
		benchmark_filtered_region=remove_IGVrm_bed.benchmark_filtered_region,
		LCL5_annotated_vcf=filter_vcf.LCL5_annotated_vcf,
		LCL6_annotated_vcf=filter_vcf.LCL6_annotated_vcf,
		LCL7_annotated_vcf=filter_vcf.LCL7_annotated_vcf,
		LCL8_annotated_vcf=filter_vcf.LCL8_annotated_vcf,
		docker=RTGdocker,
		disk_size=disk_size,
		cluster_config=cluster_config
	}

}

