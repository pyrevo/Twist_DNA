
rule all:
    input: "DATA/background.tsv"

rule Calculate_background:
    input:
        gvcfs = "DATA/gvcf_files.txt",
    output:
        background = "DATA/background.tsv",
    container:
        config["singularity"].get("python", config["singularity"].get("default", ""))
    script:
        "../../../scripts/python/Calculate_background.py"
