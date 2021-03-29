
rule msisensor2:
    input:
        bam="Bam/DNA/{sample}-ready.bam",
        bai="Bam/DNA/{sample}-ready.bam.bai",
    output:
        msi_score="Results/DNA/{sample}/MSI/{sample}.msi",
    params:
        prefix="Results/DNA/{sample}/MSI/",
        models="/opt/msisensor2/models_hg19_GRCh37",
    container:
        config["singularity"]["msisensor2"]
    shell:
        "msisensor2 msi -M {params.models} -t {input.bam} -o {params.prefix}"
