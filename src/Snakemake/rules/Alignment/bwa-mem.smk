
rule bwa_mem:

    input:
        reads=["fastq/DNA/{sample}_R1.fastq.gz", "fastq/DNA/{sample}_R2.fastq.gz"],
    output:
        bam="bam/{sample}-sort.bam",
    log:
        "logs/map/bwa/{sample}.log",
    params:
        index=config["reference"]["ref"],
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:illumina\tPU:{sample}' -v 1",
    threads: 10
    singularity:
        config["singularity"]["bwa"]
    shell:
        "(bwa mem -t {threads} {params.extra} {params.index} {input.reads} | samtools sort -o {output.bam} - ) &> {log}"


rule samtools_index:
    input:
        "bam/{sample}-sort.bam",
    output:
        "bam/{sample}-sort.bam.bai",
    log:
        "logs/map/samtools_index/{sample}.log",
    singularity:
        config["singularity"]["samtools"]
    shell:
        "(samtools index {input} {output}) &> {log}"
