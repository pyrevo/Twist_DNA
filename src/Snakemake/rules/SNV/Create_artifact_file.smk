
configfile: "Twist_DNA.yaml"


localrules:
    all,
    Create_artifacts,


rule all:
    input:
        artifacts="Artifact_positions.txt",


rule Create_artifacts:
    input:
        vcfs="DATA/Artifact_vcfs.txt",
    output:
        artifacts="Artifact_positions.txt",
    params:
        callers=["vardict", "freebayes", "mutect2", "varscan"],
    container:
        config["singularity"].get("python", config["singularity"].get("default", ""))
    script:
        "../../../scripts/python/Create_artifact_file.py"
