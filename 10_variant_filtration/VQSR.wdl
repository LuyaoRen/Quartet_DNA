version 1.0

workflow vcffilter {
    
    input {
        File vcf
        File vcf_index
        File ref_fasta = "oss://nccl-sjzp-2021/ReferenceSequence/GRCh38.d1.vd1.fa"
        Array[File] ref_fasta_index =["oss://nccl-sjzp-2021/ReferenceSequence/GRCh38.d1.vd1.dict","oss://nccl-sjzp-2021/ReferenceSequence/GRCh38.d1.vd1.fa.amb","oss://nccl-sjzp-2021/ReferenceSequence/GRCh38.d1.vd1.fa.ann","oss://nccl-sjzp-2021/ReferenceSequence/GRCh38.d1.vd1.fa.bwt","oss://nccl-sjzp-2021/ReferenceSequence/GRCh38.d1.vd1.fa.fai","oss://nccl-sjzp-2021/ReferenceSequence/GRCh38.d1.vd1.fa.pac","oss://nccl-sjzp-2021/ReferenceSequence/GRCh38.d1.vd1.fa.sa"]
        File dbsnp = "oss://pgx-reference-genome/GRCh38.d1.vd1/dbsnp_146.hg38.vcf"
        File dbsnp_index = "oss://pgx-reference-genome/GRCh38.d1.vd1/dbsnp_146.hg38.vcf.idx"
        File dbsnp_Mill = "oss://pgx-reference-genome/GRCh38.d1.vd1/Mills_and_1000G_gold_standard.indels.hg38.vcf"
        File dbsnp_Mill_index = "oss://pgx-reference-genome/GRCh38.d1.vd1/Mills_and_1000G_gold_standard.indels.hg38.vcf.idx"
        File dbsnp_1000G_omni = "oss://pgx-reference-genome/GRCh38.d1.vd1/1000G_omni2.5.hg38.vcf.gz"
        File dbsnp_1000G_omni_index = "oss://pgx-reference-genome/GRCh38.d1.vd1/1000G_omni2.5.hg38.vcf.gz.tbi"
        File dbsnp_hapmap = "oss://pgx-reference-genome/GRCh38.d1.vd1/hapmap_3.3.hg38.vcf.gz"
        File dbsnp_hapmap_index = "oss://pgx-reference-genome/GRCh38.d1.vd1/hapmap_3.3.hg38.vcf.gz.tbi"
        File dbsnp_1000G_phase1 = "oss://pgx-reference-genome/GRCh38.d1.vd1/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
        File dbsnp_1000G_phase1_index = "oss://pgx-reference-genome/GRCh38.d1.vd1/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi"
        String sample
    }
    
    # 调用任务或者子流程, 提供对应的输入参数
    # https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#call-statement
    call vqsr {
        input:
            vcf = vcf,
            vcf_index = vcf_index,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            dbsnp = dbsnp,
            dbsnp_index = dbsnp_index,
            dbsnp_Mill = dbsnp_Mill,
            dbsnp_Mill_index = dbsnp_Mill_index,
            dbsnp_1000G_omni = dbsnp_1000G_omni,
            dbsnp_1000G_omni_index = dbsnp_1000G_omni_index,
            dbsnp_hapmap = dbsnp_hapmap,
            dbsnp_hapmap_index = dbsnp_hapmap_index,
            dbsnp_1000G_phase1 = dbsnp_1000G_phase1,
            dbsnp_1000G_phase1_index = dbsnp_1000G_phase1_index,
            sample=sample
    }
    
    # 流程分析最终结果
    # https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#outputs
    output {
        File masked_vcf = vqsr.masked_vcf
        File filtered_vcf = vqsr.filtered_vcf
    }
}

# 在此定义分析流程中的任务
# https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#task-definition

task vqsr {

    # 必选，任务输入参数，支持多种类型变量声明。
    # https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#types

    input {
        File vcf
        File vcf_index
        File ref_fasta
        Array[File] ref_fasta_index
        File dbsnp
        File dbsnp_index
        File dbsnp_Mill
        File dbsnp_Mill_index
        File dbsnp_1000G_omni
        File dbsnp_1000G_omni_index
        File dbsnp_hapmap
        File dbsnp_hapmap_index
        File dbsnp_1000G_phase1
        File dbsnp_1000G_phase1_index
        String sample

        Int cpu = 32
        String memory = "64G"
        String disks = "local-disk 250 cloud_ssd"
    }



    # 任务命令行，为工具的实际执行脚本
    # https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#command-section

    command <<<
        resource_text="--resource ~{dbsnp_1000G_phase1} \
              --resource_param 1000G,known=false,training=true,truth=false,prior=10.0 "
        resource_text="$resource_text --resource ~{dbsnp_1000G_omni} \
              --resource_param omni,known=false,training=true,truth=true,prior=12.0 "
        resource_text="$resource_text --resource ~{dbsnp} \
              --resource_param dbsnp,known=true,training=false,truth=false,prior=2.0 "
        resource_text="$resource_text --resource ~{dbsnp_hapmap} \
              --resource_param hapmap,known=false,training=true,truth=true,prior=15.0"

        annotation_array="QD DP FS SOR MQ ReadPosRankSum MQRankSum"
        annotate_text=""
        for annotation in $annotation_array; do
          annotate_text="$annotate_text --annotation $annotation"
        done

        sentieon driver -r ~{ref_fasta} --algo VarCal -v ~{vcf} $resource_text $annotate_text --var_type SNP --plot_file ~{sample}.vqsrSNP.hc.plotfile --tranches_file ~{sample}.vqsrSNP.hc.tranches ~{sample}.vqsrSNP.hc.recal

        sentieon driver -r ~{ref_fasta} --algo ApplyVarCal -v ~{vcf} --var_type SNP --tranches_file ~{sample}.vqsrSNP.hc.tranches --sensitivity 99.0 --recal ~{sample}.vqsrSNP.hc.recal ~{sample}.vqsrSNP.hc.recaled.vcf.gz

        sentieon plot vqsr -o ~{sample}.vqsrSNP.pdf ~{sample}.vqsrSNP.hc.plotfile

        resource_text="--resource ~{dbsnp_Mill} \
              --resource_param Mills,known=false,training=true,truth=true,prior=12.0 "
        resource_text="$resource_text --resource ~{dbsnp} \
              --resource_param dbsnp,known=true,training=false,truth=false,prior=2.0 "

        annotation_array="QD DP FS SOR ReadPosRankSum MQRankSum"
        annotate_text=""
        for annotation in $annotation_array; do
          annotate_text="$annotate_text --annotation $annotation"
        done

        sentieon driver -r ~{ref_fasta} --algo VarCal -v ~{vcf} $resource_text $annotate_text --var_type INDEL --plot_file ~{sample}.vqsrINDEL.hc.plotfile --max_gaussians 4 --tranches_file ~{sample}.vqsrINDEL.hc.tranches ~{sample}.vqsrINDEL.hc.recal

        sentieon driver -r ~{ref_fasta} --algo ApplyVarCal -v ~{sample}.vqsrSNP.hc.recaled.vcf.gz --var_type INDEL --recal ~{sample}.vqsrINDEL.hc.recal --tranches_file ~{sample}.vqsrINDEL.hc.tranches --sensitivity 99.0 ~{sample}.vqsrSNPINDEL.hc.recaled.vcf.gz

        sentieon plot vqsr -o ~{sample}.vqsrINDEL.VQSR.pdf ~{sample}.vqsrINDEL.hc.plotfile

        zcat ~{sample}.vqsrSNPINDEL.hc.recaled.vcf.gz | grep '#' > header
        zcat ~{sample}.vqsrSNPINDEL.hc.recaled.vcf.gz | grep -v '#' | grep PASS > body
        cat header body > ~{sample}.vqsrSNPINDEL.hc.recaled.pass.vcf

    >>>
    
    # 可选，默认提供1核2G，40G硬盘的执行环境
    # https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#runtime-section
    # https://help.aliyun.com/document_detail/257724.html
    runtime {
        cpu: cpu
        memory: memory
        disks: disks
        software: "sentieon:202112.05"
    }

    # 任务运行输出结果
    # https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#outputs-section
    output {
        File masked_vcf = "~{sample}.vqsrSNPINDEL.hc.recaled.vcf.gz"
        File filtered_vcf = "~{sample}.vqsrSNPINDEL.hc.recaled.pass.vcf"
    }
    
}
