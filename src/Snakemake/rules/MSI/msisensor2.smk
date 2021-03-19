
rule msisensor2:
    input:
        bam="Bam/DNA/{sample}-ready.bam",
        bai="Bam/DNA/{sample}-ready.bam.bai",
    output:
        msi_score="Results/DNA/{sample}/MSI/{sample}.msi",
    params:
        prefix="MSI/{sample}.msi",
        models="/opt/msisensor2/models_hg38",
    singularity:
        config["singularity"]["msisensor2"]
    shell:
        "msisensor2 msi -M {params.models} -t {input.bam} -o {params.prefix}"
