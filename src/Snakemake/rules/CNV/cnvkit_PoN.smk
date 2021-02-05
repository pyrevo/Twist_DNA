
configfile: "Twist_DNA.yaml"


rule all:
    input:
        PoN="DATA/cnvkit_Twist_PoN.cnn",


rule Create_targets:
    input:
        bed=config["bed"]["bedfile"],
    output:
        bed="CNV/bed/cnvkit_manifest.target.bed",
    log:
        "logs/CNV_cnvkit/Create_targets.log",
    singularity:
        config["singularity"]["cnvkit"]
    shell:
        "(cnvkit.py target --split {input.bed} -o {output.bed}) &> {log}"


rule Create_anti_targets:
    input:
        bed="CNV/bed/cnvkit_manifest.target.bed",
    output:
        bed="CNV/bed/cnvkit_manifest.antitarget.bed",
    log:
        "logs/CNV_cnvkit/Create_anti_targets.log",
    singularity:
        config["singularity"]["cnvkit"]
    shell:
        "(cnvkit.py antitarget {input.bed} -o {output.bed}) &> {log}"


rule Build_normal_reference:
    input:
        #bams=expand("{normal_sample}", normal_sample=config["DNA_Samples"]),
        bams=["DNA_bam/" + s + "-ready.bam" for s in config["DNA_Samples"]],
        bed1="CNV/bed/cnvkit_manifest.target.bed",
        bed2="CNV/bed/cnvkit_manifest.antitarget.bed",
        ref=config["reference"]["ref"],
        mappability="DATA/access-5k-mappable.hg19.bed",
    output:
        PoN="DATA/cnvkit_Twist_PoN.cnn",
    threads: 4
    singularity:
        config["singularity"]["cnvkit"]
    shell:
        "cnvkit.py batch "
        "-n {input.bams} "
        "-m hybrid "
        "--output-reference {output.PoN} "
        "-t {input.bed1} "
        "-f {input.ref} "
        "-a {input.bed2} "
        "-g {input.mappability} "
        "-p {threads}"
