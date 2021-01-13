

rule GPU_align:
    input:
        fastq1 = "fastq/DNA/{sample}_R1.fastq.gz"
        fastq2 = "fastq/DNA/{sample}_R2.fastq.gz"
    output:
        bam = "DNA_bam/{sample}-ready.bam"
    params:
        ref = config["reference"]["ref"]
    log:
        "logs/GPU_align/{sample}.log"
    threads:
        40
    singularity:
        config["singularitys"]["GPU"]
    shell:
        "(pbrun fq2bam "
        "--ref {params.ref} "
        "--in-fq {input.fastq1} {input.fastq2} "
        "--out-bam {output.bam} "
        "-tmp-dir GPU_run "
        ") &> {log}"
