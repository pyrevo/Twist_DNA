
rule gridss:
    input:
        bam="Bam/DNA/{sample}-ready.bam",
        bai="Bam/DNA/{sample}-ready.bam.bai",
    output:
        vcf="Results/DNA/{sample}/gridss/{sample}.vcf.gz",
        bam="Results/DNA/{sample}/gridss/{sample}.bam",
    params:
        ref=config["reference"]["ref"],
        blacklist="DATA/hg19-blacklist.v2.bed",
        workingdir="gridss/wd/",
    threads: 5
    log:
        "logs/DNA_fusion/gridss/{sample}.log",
    singularity:
        config["singularity"].get("gridss", config["singularity"].get("default", ""))
    shell:
        "(gridss.sh --reference {params.ref} --output {output.vcf} --assembly {output.bam} --threads {threads} --blacklist {params.blacklist} --workingdir {params.workingdir} {input.bam}) > {log}"
