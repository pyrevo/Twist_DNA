

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
        OutputPath="Results/DNA/{sample}/JuLI",
        sample_name=lambda wildcards: wildcards.sample,
    threads: 10
    log:
        "logs/DNA_fusion_call/JuLI/{sample}.log",
    container:
        config["singularity"].get("JuLI", config["singularity"].get("default", ""))
    shell:
        "(Rscript -e '"
        "library(juliv0.1.6.1); "
        "callfusion(CaseBam=\"{input.bam}\", "
        "TestID=\"{params.sample_name}\", "
        "OutputPath=\"{params.OutputPath}\", "
        "Thread=\"{threads}\", "
        "Refgene=\"{params.Refgene}\", "
        "Gap=\"{params.Gap}\", "
        "Reference=\"{params.ref}\")')  &> {log}"


rule JuLI_annotate:
    input:
        fusions="Results/DNA/{sample}/JuLI/{sample}.txt",
    output:
        fusions="Results/DNA/{sample}/JuLI/{sample}.annotated.txt",
    params:
        ref=config["reference"]["ref"],
        Refgene="/opt/references/refGene_hg19.txt",
        Cosmic="/opt/references/CosmicFusionExport_V76.tsv",
        Pfam="/opt/references/Pfam-A.full.human",
        Uniprot="/opt/references/HGNC_GeneName_UniProtID_160524.txt",
    log:
        "logs/DNA_fusion_annotate/JuLI/{sample}.log",
    container:
        config["singularity"].get("JuLI", config["singularity"].get("default", ""))
    shell:
        "(Rscript -e '"
        "library(juliv0.1.6.1); "
        "annofusion(Output=\"{input.fusions}\", "
        "Refgene=\"{params.Refgene}\", "
        "Cosmic=\"{params.Cosmic}\", "
        "Pfam=\"{params.Pfam}\", "
        "Uniprot=\"{params.Uniprot}\")')  &> {log}"
