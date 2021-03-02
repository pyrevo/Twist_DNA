rule samtools_stats:
    input:
        bam="DNA_bam/{sample}-ready.bam",
    output:
        "qc/{sample}/{sample}.samtools-stats.txt",
    params:
        extra="-t " + config["bed"]["bedfile"],
    log:
        "logs/qc/samtools_stats/{sample}.log",
    singularity:
        config["singularity"]["samtools"]
    shell:
        "(samtools stats {params.extra} {input} > {output} ) &> {log}"
