
rule fixAF:
    input:
        "{method}/{sample}.{method}.fixAF.vcf",
    output:
        temp("{method}/{sample}.{method}.okAF.vcf"),
    log:
        "logs/variantCalling/fixAF/{method}/{sample}.log",
    container:
        config["singularity"].get("python", config["singularity"].get("default", ""))
    script:
        "../../../scripts/python/fix_af.py"
