
rule move_fastq:
    input:
        fastq1="fastq_temp/DNA/{sample}_R1.fastq.gz",
        fastq2="fastq_temp/DNA/{sample}_R2.fastq.gz",
    output:
        fastq1="fastq/DNA/{sample}_R1.fastq.gz",
        fastq2="fastq/DNA/{sample}_R2.fastq.gz",
    shell:
        "mv {input.fastq1} {output.fastq1} && "
        "mv {input.fastq2} {output.fastq2}"
