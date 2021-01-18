
configfile: "Twist_DNA.yaml"

rule all:
    input:
        PoN="DATA/cnvkit_PoN.cnn",

rule Build_normal_reference:
    input:
        bams = expand("{normal_sample}", normal_sample=config["Normal_samples"]),
        #bams=["DNA_bam/" + s + "-ready.bam" for s in config["DNA_Samples"]],
        bed1="bed/manifest.target.bed",
        bed2="bed/manifest.antitarget.bed",
        ref=config["reference"]["ref"],
        mappability = "DATA/access-5k-mappable.hg19.bed",
    output:
        PoN="DATA/cnvkit_PoN.cnn",
    threads: 4
    singularity: config["singularity"]["cnvkit"]
    shell:
        "cnvkit.py batch "
        "-n {input.bams} "
        "-m hybrid "
        "--output-reference {output.PoN} "
        "-t {input.bed1} "
        "-f {input.reference} "
        "-a {input.bed2} "
        "-g {input.mappability} "
        "-p {threads}"
