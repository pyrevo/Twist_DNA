
rule geneFuse:
    input:
        fastq1="fastq/DNA/{sample}_R1.fastq.gz",
        fastq2="fastq/DNA/{sample}_R2.fastq.gz",
    output:
        html="Results/DNA/{sample}/geneFuse/geneFuse_report_{sample}.html",
        fusions="Results/DNA/{sample}/geneFuse/fusions_{sample}.txt",
    params:
        genes=config["geneFuse"]["genes"],
        ref=config["reference"]["ref"],
    threads: 4
    log:
        "logs/DNA_fusion/geneFuse/{sample}.log",
    container:
        config["singularity"].get("geneFuse", config["singularity"].get("default", ""))
    shell:
        "(genefuse -r {params.ref} -f {params.genes} -1 {input.fastq1} -2 {input.fastq2} -h {output.html} > {output.fusions}) &> {log}"
