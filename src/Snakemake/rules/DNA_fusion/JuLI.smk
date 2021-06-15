

rule JuLI_call:
    input:
        bam="Bam/DNA/{sample}-ready.bam",
        bai="Bam/DNA/{sample}-ready.bam.bai",
    output:
        fusions="Results/DNA/{sample}/JuLI/{sample}.txt",
    params:
        ref=config["reference"]["ref"],
        Refgene="/opt/references/refGene_hg19.txt",
        Gap="/opt/references/gap_hg19.txt",
        OutputPath="Results/DNA/{sample}/JuLI/{sample}",
    threads: 10
    log:
        "logs/DNA_fusion/JuLI/{sample}.log",
    container:
        config["singularity"].get("JuLI", config["singularity"].get("default", ""))
    shell:
        "Rscript -e '"
        "callfusion(CaseBam=\"{input.bam}\", "
        "TestID='{sample}', "
        "OutputPath='{sample}', "
        "Thread={threads}, "
        "Refgene=\"{params.Refgene}\", "
        "Gap=\"{params.Gap}\", "
        "Reference=\"{params.ref}\")'"


rule JuLI_annotate:
    input:
        fusions="Results/DNA/{sample}/JuLI/{sample}.txt",
    output:
        fusions="Results/DNA/{sample}/JuLI/{sample}.annotation.txt",
    params:
        ref=config["reference"]["ref"],
        Refgene="/opt/references/refGene_hg19.txt",
        Cosmic="/opt/references/CosmicFusionExport_V76.tsv",
        Pfam="/opt/references/Pfam-A.full.human",
        Uniprot="/opt/references/HGNC_GeneName_UniProtID_160524.txt",
    log:
        "logs/DNA_fusion/JuLI/{sample}.log",
    container:
        config["singularity"].get("JuLI", config["singularity"].get("default", ""))
    shell:
        "Rscript -e '"
        "annofusion(Output=\"{input.fusions}\", "
        "Refgene=\"{params.Refgene}\", "
        "Cosmic=\"{params.Cosmic}\", "
        "Pfam=\"{params.Pfam}\", "
        "Uniprot=\"{params.Uniprot}\")"
