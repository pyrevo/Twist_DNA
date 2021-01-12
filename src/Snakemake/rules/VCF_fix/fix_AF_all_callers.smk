
rule fixAF:
    input:
        "{method}/{sample}.{method}.fixAF.vcf"
    output:
        temp("{method}/{sample}.{method}.okAF.vcf")
    log:
        "logs/variantCalling/fixAF/{method}/{sample}.log"
    singularity:
        config["singularity"]["python"]
    shell:
        "(python3.6 src/Snakemake/scripts/fix_af.py {input} {output}) &> {log}"

localrules: bgzipCallers
rule bgzipCallers:
    input:
        vcf = "{method}/{sample}.{method}.okAF.vcf"
    output:
        vcf = "{method}/{sample}.{method}.okAF.vcf.gz",
        tabix = "{method}/{sample}.{method}.okAF.vcf.gz.tbi"
    log:
        "logs/variantCalling/bgzip/{method}/{sample}.log"
    singularity:
        config["singularity"]["bcftools"]
    shell:
        "(bgzip {input.vcf} && tabix {output.vcf}) 2> {log}"
