
rule Find_multibp_SNV:
    input:
        vcf="recall/{sample}.ensemble.vep.exon.soft_filter.vcf",
        ref=config["reference"]["ref"],
    output:
        vcf=temp("recall/{sample}.ensemble.vep.exon.soft_filter.multibp.vcf.temp"),
    singularity:
        config["singularity"]["python"]
    script:
        "../../../scripts/python/Multibp_SNV.py"

rule sort_multiplebp_vcf:
    input:
        vcf="recall/{sample}.ensemble.vep.exon.soft_filter.multibp.vcf.temp",
    output:
        vcf="recall/{sample}.ensemble.vep.exon.soft_filter.multibp.vcf",
    singularity:
        config["singularity"]["bcftools"]
    shell:
        "bcftools sort -o {output.vcf} -O v {input.vcf}"
