
configfile: "Twist_DNA.yaml"

localrules:
    all,
    Create_varscan_artifacts,

rule all:
    input:
        artifacts="Artifact_positions.txt",


rule Create_varscan_artifacts:
    input:
        # {path}/{sample}.ensemble.vep.exon.soft_filter.multibp.vcf
        vcfs="Artifact_vcfs.txt",
    output:
        artifacts="Artifact_positions.txt",
    params:
        callers=["Vardict", "Freebayes", "Mutect2", "Varscan"],
    container:
        config["singularity"].get("python", config["singularity"].get("default", ""))
    script:
        "../../../scripts/python/Create_artifact_file.py"
