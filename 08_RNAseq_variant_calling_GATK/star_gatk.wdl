version 1.0

# 在此定义您的分析流程
workflow RNAseq {
    
    # 流程输入文件和参数
    input {
        File fastq1
        File fastq2
        String sample_id
        String platform = "Illumina"
        File ref_fasta = "oss://nccl-sjzp-2021/ReferenceSequence/GRCh38.d1.vd1.fa"
        Array[File] ref_fasta_index =["oss://nccl-sjzp-2021/ReferenceSequence/GRCh38.d1.vd1.dict","oss://nccl-sjzp-2021/ReferenceSequence/GRCh38.d1.vd1.fa.amb","oss://nccl-sjzp-2021/ReferenceSequence/GRCh38.d1.vd1.fa.ann","oss://nccl-sjzp-2021/ReferenceSequence/GRCh38.d1.vd1.fa.bwt","oss://nccl-sjzp-2021/ReferenceSequence/GRCh38.d1.vd1.fa.fai","oss://nccl-sjzp-2021/ReferenceSequence/GRCh38.d1.vd1.fa.pac","oss://nccl-sjzp-2021/ReferenceSequence/GRCh38.d1.vd1.fa.sa"]
        File dbsnp_vcf = "oss://pgx-reference-genome/GRCh38.d1.vd1/dbsnp_146.hg38.vcf"
        File dbsnp_index = "oss://pgx-reference-genome/GRCh38.d1.vd1/dbsnp_146.hg38.vcf.idx"
        File dbmills_vcf = "oss://pgx-reference-genome/GRCh38.d1.vd1/Mills_and_1000G_gold_standard.indels.hg38.vcf"
        File dbmills_index = "oss://pgx-reference-genome/GRCh38.d1.vd1/Mills_and_1000G_gold_standard.indels.hg38.vcf.idx"
        Array[File] star_index = ["oss://pgx-reference-genome/STAR_index/Genome","oss://pgx-reference-genome/STAR_index/SA","oss://pgx-reference-genome/STAR_index/SAindex","oss://pgx-reference-genome/STAR_index/chrLength.txt","oss://pgx-reference-genome/STAR_index/chrName.txt","oss://pgx-reference-genome/STAR_index/chrNameLength.txt","oss://pgx-reference-genome/STAR_index/chrStart.txt","oss://pgx-reference-genome/STAR_index/exonGeTrInfo.tab","oss://pgx-reference-genome/STAR_index/exonInfo.tab","oss://pgx-reference-genome/STAR_index/geneInfo.tab","oss://pgx-reference-genome/STAR_index/genomeParameters.txt","oss://pgx-reference-genome/STAR_index/sjdbInfo.txt","oss://pgx-reference-genome/STAR_index/sjdbList.fromGTF.out.tab","oss://pgx-reference-genome/STAR_index/sjdbList.out.tab","oss://pgx-reference-genome/STAR_index/transcriptInfo.tab"]
        File bed = "oss://pgx-dna-results/reference_dataset/Quartet.high.confidence.region.v202103.removedSV.cds.bed"
    }
    
    # 调用子任务

    call SentieonSTARToVcf {
        input:
            fastq1=fastq1,
            fastq2=fastq2,
            ref_fasta=ref_fasta,
            ref_fasta_index=ref_fasta_index,
            dbsnp_vcf=dbsnp_vcf,
            dbsnp_index=dbsnp_index,
            dbmills_vcf=dbmills_vcf,
            dbmills_index=dbmills_index,
            sample_id=sample_id,
            star_index=star_index,
            platform=platform,
            bed=bed
    }
    
    # 流程分析结果
    output {
        Pair[File,File] bam = SentieonSTARToVcf.bam
        Pair[File,File] vcf = SentieonSTARToVcf.vcf
        Array[File] star_out = SentieonSTARToVcf.star_out
    }
}

# 在此定义流程中使用的工具
task SentieonSTARToVcf {
    # 工具输入文件和参数
    input {
        File fastq1
        File fastq2
        String sample_id
        String platform
        File ref_fasta
        Array[File] ref_fasta_index
        File dbsnp_vcf
        File dbsnp_index
        File dbmills_vcf
        File dbmills_index
        Array[File] star_index
        File bed

        Int cpu = 32
        String memory = "64G"
        String disks = "local-disk 250 cloud_ssd"
    }
    
    String out_bam = sample_id + ".splitted.bam"
    String out_bai = sample_id + ".splitted.bam.bai"
    String out_vcf = sample_id + "_hc.vcf"
    String out_vcf_idx = sample_id + "_hc.vcf.idx"

    # 工具运行命令
    command <<<
        set -exo pipefail

        mkdir genomeDir
        for i in ~{sep="\t" star_index}; do cp $i genomeDir; done
        
        sentieon STAR \
            --twopassMode Basic --twopass1readsN -1 --genomeDir genomeDir \
            --runThreadN ~{cpu} \
            --outStd BAM_Unsorted --outSAMtype BAM Unsorted \
            --outBAMcompression 0 \
            --readFilesIn ~{fastq1} ~{fastq2} \
            --readFilesCommand "zcat" \
            --outFileNamePrefix STARout \
            --sjdbOverhang 149 \
            --outSAMattrRGline ID:~{sample_id} SM:~{sample_id} PL:~{platform} | \
        sentieon util sort -i - -r ~{ref_fasta} -t ~{cpu} -o ~{sample_id}.sorted.bam

        sentieon driver -r ~{ref_fasta} -t ~{cpu} -i ~{sample_id}.sorted.bam --traverse_param 1000000/10000 --algo LocusCollector --fun score_info ~{sample_id}.score.txt.gz

        sentieon driver -r ~{ref_fasta} -t ~{cpu} -i ~{sample_id}.sorted.bam --traverse_param 1000000/10000 \
         --algo Dedup --rmdup --score_info ~{sample_id}.score.txt.gz --metrics ~{sample_id}.dedup_metrics.txt ~{sample_id}.deduped.bam

        sentieon driver -r ~{ref_fasta} -t ~{cpu} -i ~{sample_id}.deduped.bam --algo RNASplitReadsAtJunction --reassign_mapq 255:60 ~{out_bam}

        sentieon driver --traverse_param 1000000/10000 -r ~{ref_fasta} -t ~{cpu} --interval ~{bed} -i ~{out_bam} --algo QualCal -k ~{dbsnp_vcf} -k ~{dbmills_vcf} ~{sample_id}_recal_data.table

        sentieon driver --traverse_param 1000000/10000 -r ~{ref_fasta} -t ~{cpu} -i ~{out_bam} --interval ~{bed} -q ~{sample_id}_recal_data.table --algo QualCal -k ~{dbsnp_vcf} -k ~{dbmills_vcf} ~{sample_id}_recal_data.table.post --algo ReadWriter ~{sample_id}.sorted.deduped.recaled.bam

        sentieon driver --traverse_param 1000000/10000 -t ~{cpu} --interval ~{bed} --algo QualCal --plot --before ~{sample_id}_recal_data.table --after ~{sample_id}_recal_data.table.post ~{sample_id}_recal_data.csv

        sentieon driver --traverse_param 1000000/10000 -r ~{ref_fasta} --interval ~{bed} -t ~{cpu} -i ~{sample_id}.sorted.deduped.recaled.bam --algo Haplotyper -d ~{dbsnp_vcf} --trim_soft_clip --call_conf 20 --emit_conf 20 ~{sample_id}_hc.vcf

    >>>
    runtime {
        cpu: cpu
        memory: memory
        disks: disks
        software: "sentieon:202112.05"
    }
    # 工具运行输出结果
    output {
        Array[File] star_out = glob("STARout*")
        Pair[File,File] bam = (out_bam,out_bai)
        Pair[File,File] vcf = (out_vcf,out_vcf_idx)
    }
}