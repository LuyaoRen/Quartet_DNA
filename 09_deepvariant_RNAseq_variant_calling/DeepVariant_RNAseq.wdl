version 1.0

workflow RNAvariants_DV {
    
    input {
        File bam
        File bai
        File bed
        File ref_fasta = "oss://nccl-sjzp-2021/ReferenceSequence/GRCh38.d1.vd1.fa"
        Array[File] ref_fasta_index =["oss://nccl-sjzp-2021/ReferenceSequence/GRCh38.d1.vd1.dict","oss://nccl-sjzp-2021/ReferenceSequence/GRCh38.d1.vd1.fa.amb","oss://nccl-sjzp-2021/ReferenceSequence/GRCh38.d1.vd1.fa.ann","oss://nccl-sjzp-2021/ReferenceSequence/GRCh38.d1.vd1.fa.bwt","oss://nccl-sjzp-2021/ReferenceSequence/GRCh38.d1.vd1.fa.fai","oss://nccl-sjzp-2021/ReferenceSequence/GRCh38.d1.vd1.fa.pac","oss://nccl-sjzp-2021/ReferenceSequence/GRCh38.d1.vd1.fa.sa"]
        String sample
        Array[File] models = ["oss://pgx-reference-genome/deepvariant_model/rna_seq/models_DeepVariant_1.4.0_DeepVariant-inception_v3-1.4.0+data-rnaseq_standard_model.ckpt.data-00000-of-00001","oss://pgx-reference-genome/deepvariant_model/rna_seq/models_DeepVariant_1.4.0_DeepVariant-inception_v3-1.4.0+data-rnaseq_standard_model.ckpt.example_info.json","oss://pgx-reference-genome/deepvariant_model/rna_seq/models_DeepVariant_1.4.0_DeepVariant-inception_v3-1.4.0+data-rnaseq_standard_model.ckpt.index","oss://pgx-reference-genome/deepvariant_model/rna_seq/models_DeepVariant_1.4.0_DeepVariant-inception_v3-1.4.0+data-rnaseq_standard_model.ckpt.meta"]
    }
    
    # 调用任务或者子流程, 提供对应的输入参数
    # https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#call-statement
    call deepvariant {
        input:
            bam = bam,
            bed = bed,
            bai = bai,
            models = models,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            sample=sample
    }
    
    output {
        File vcf = deepvariant.vcf
    }
}

task deepvariant {

    input {
        File bam
        File bed
        File bai
        File ref_fasta
        Array[File] ref_fasta_index
        Array[File] models
        String sample

        Int cpu = 32
        String memory = "64G"
        String disks = "local-disk 250 cloud_ssd"
    }

    command <<<
        mkdir model
        for i in ~{sep="\t" models}; do cp $i model; done

        /opt/deepvariant/bin/run_deepvariant \
            --model_type=WES \
            --customized_model=model/models_DeepVariant_1.4.0_DeepVariant-inception_v3-1.4.0+data-rnaseq_standard_model.ckpt \
            --ref=~{ref_fasta} \
            --reads=~{bam} \
            --output_vcf=~{sample}.deepvariant.vcf.gz \
            --num_shards=~{cpu} \
            --regions=~{bed}\
            --make_examples_extra_args="split_skip_reads=true,channels=''" \
            --intermediate_results_dir output/intermediate_results_dir

    >>>

    runtime {
        cpu: cpu
        memory: memory
        disks: disks
        docker: "registry.cn-shanghai.aliyuncs.com/pgx-dna/deepvariant:202306"
    }
    output {
        File vcf = "~{sample}.deepvariant.vcf.gz"
    }
    
}
