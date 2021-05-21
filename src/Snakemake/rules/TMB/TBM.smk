
rule TMB:
    input:
        vcf="Results/DNA/{sample}/vcf/{sample}.ensemble.vep.exon.soft_filter.multibp.vcf",
        artifacts="DATA/TMB_artifact_positions.txt",
    output:
        tmb="Results/DNA/{sample}/TMB/{sample}.TMB.txt",
    container:
        config["singularity"].get("python", config["singularity"].get("default", ""))
    script:
        "../../../scripts/python/TMB.py"
