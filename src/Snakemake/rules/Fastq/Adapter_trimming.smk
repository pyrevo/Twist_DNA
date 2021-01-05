

rule cutadapt:
    input:
        fastq1 = "fastq_temp/DNA/{sample}_R1.fastq.gz",
        fastq2 = "fastq_temp/DNA/{sample}_R2.fastq.gz"
    output:
        fastq1 = "fastq/DNA/{sample}_R1.fastq.gz",
        fastq2 = "fastq/DNA/{sample}_R2.fastq.gz",
        qc = "fastq/DNA/{sample}.qc.txt"
    params:
         # https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types
        adapters_r1 = "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",  #"-a AGAGCACACGTCTGAACTCCAGTCAC -g AGATCGGAAGAGCACACGT",
        adapters_r2 = "-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT", #"-A AGAGCACACGTCTGAACTCCAGTCAC -G AGATCGGAAGAGCACACGT",
        # https://cutadapt.readthedocs.io/en/stable/guide.html#
        others = "--minimum-length 2 -q 20"
    log:
        "logs/trimming/cutadapt/{sample}.log"
    threads: 10
    singularity:
        config["singularity"]["cutadapt"]
    shell:
        "(cutadapt {params.adapters_r1} {params.adapters_r2} {params.others} -o {output.fastq1} -p {output.fastq2} -j {threads} {input.fastq1} {input.fastq2} > {output.qc} ) &> {log}"
