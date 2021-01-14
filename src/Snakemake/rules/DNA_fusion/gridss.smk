
rule gridss:
    input:
        bam = "DNA_bam/{sample}-ready.bam",
        bai = "DNA_bam/{sample}-ready.bam.bai"
    output:
        vcf = "gridss/{sample}.vcf.gz",
        bam = "gridss/{sample}.bam"
    params:
        ref = config["reference"]["ref"],
        blacklist = "Data/hg19-blacklist.v2.bed"
    threads: 5
    log:
        "logs/DNA_fusion/gridss/{sample}.log"
    singularity:
        config["singularity"]["gridss"]
    shell:
        "(gridss.sh --reference {params.ref} --output {output.vcf} --assembly {output.bam} --threads {threads} --blacklist {params.blacklist} {input.bam}) > {log}"
