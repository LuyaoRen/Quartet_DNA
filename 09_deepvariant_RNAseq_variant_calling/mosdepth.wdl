version 1.0

workflow BAMcov {
    
    input {
        File bam
        File bai
        String sample
    }
    
    call mosdepth {
        input:
            bam = bam,
            bai = bai,
            sample=sample
    }
    
    output {
        Array[File] coverage = mosdepth.coverage
    }
}

task mosdepth {

    input {
        File bam
        File bai
        String sample

        Int cpu = 32
        String memory = "64G"
        String disks = "local-disk 250 cloud_ssd"
    }

    command <<<

        mosdepth --threads ~{cpu} ~{sample}.mosdepth ~{bam}

    >>>

    runtime {
        cpu: cpu
        memory: memory
        disks: disks
        docker:"registry.cn-shanghai.aliyuncs.com/pgx-dna/mosdepth:0.3.3"
    }
    output {
        Array[File] coverage = glob("*.mosdepth*")
    }
    
}
