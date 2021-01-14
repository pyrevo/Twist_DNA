


rule fastp:
    input:
        fastq1 = "fastq_temp/DNA/{sample}_R1.fastq.gz",
        fastq2 = "fastq_temp/DNA/{sample}_R2.fastq.gz"
    output:
        fastq1 = "fastq/DNA/{sample}_R1.fastq.gz",
        fastq2 = "fastq/DNA/{sample}_R2.fastq.gz"
    params:
        html = "fastq/DNA/{sample}.html",
        json = "fastq/DNA/{sample}.json"
    log:
        "logs/trimming/fastp/{sample}.log"
    threads: 5
    singularity:
        config["singularity"]["fastp"]
    shell:
        "(fastp -i {input.fastq1} -I {input.fastq2} -o {output.fastq1} -O {output.fastq2} -w {threads} -h {params.html} -j {params.json}) &> {log}"
