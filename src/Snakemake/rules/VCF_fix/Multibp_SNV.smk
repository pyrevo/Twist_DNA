
rule Find_multibp_SNV:
    input:
        vcf="recall/{sample}.ensemble.final.exon.soft_filter.ffpe.vcf",
        ref=config["reference"]["ref"],
    output:
        vcf="recall/{sample}.ensemble.final.exon.soft_filter.ffpe.multibp.vcf",
    singularity:
        config["singularity"]["python"]
    script:
        "../../../scripts/python/Multibp_SNV.py"
