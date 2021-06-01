
rule msisensor_pro:
    input:
        bam="Bam/DNA/{sample}-ready.bam",
        bai="Bam/DNA/{sample}-ready.bam.bai",
        PoN=config["PoN"]["msisensor-pro"],
    output:
        msi_score="/MSI/{sample}",
    params:
        extra="-c 50 -b 2",  # -c = minimal coverage
        out_prefix="/MSI/{sample}",
    log:
        "logs/MSI/msisensor-pro_{sample}.log",
    container:
        config["singularity"].get("msisensor-pro", config["singularity"].get("default", ""))
    shell:
        "(msisensor-pro pro {params.extra} -d {input.PoN} -t {input.bam} -o {params.out_prefix}) > {log}"


rule msisensor_pro_copy_results:
    input:
        msi_score="/MSI/{sample}",
    output:
        msi_score="Results/DNA/{sample}/MSI/{sample}",
    log:
        "logs/MSI/msisensor-pro_copy_{sample}.log",
    shell:
        "(cp {input.msi_score} {output.msi_score}) > {log}"
