
rule Create_targets:
    input:
        bed=config["bed"]["bedfile"],
    output:
        bed="CNV/bed/cnvkit_manifest.target.bed",
    log:
        "logs/CNV_cnvkit/Create_targets.log",
    singularity:
        config["singularity"].get("cnvkit", config["singularity"].get("default", ""))
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
        config["singularity"].get("cnvkit", config["singularity"].get("default", ""))
    shell:
        "(cnvkit.py antitarget {input.bed} -o {output.bed}) &> {log}"


rule Call_cnv:
    input:
        bams=["Bam/DNA/" + s + "-ready.bam" for s in config["DNA_Samples"]],
        PoN=config["PoN"]["cnvkit"],
    output:
        regions=["CNV/cnvkit_calls/" + sample_id + "-ready.cnr" for sample_id in config["DNA_Samples"]],
        segments=["CNV/cnvkit_calls/" + sample_id + "-ready.cns" for sample_id in config["DNA_Samples"]],
    params:
        outdir="CNV/cnvkit_calls/",
    log:
        "logs/CNV_cnvkit/Call_cnv.log",
    threads: 8
    singularity:
        config["singularity"].get("cnvkit", config["singularity"].get("default", ""))
    shell:
        "(cnvkit.py batch {input.bams} -r {input.PoN} -p {threads} -d {params.outdir}) &> {log}"


rule Filter_cnv:
    input:
        segments=["CNV/cnvkit_calls/" + sample_id + "-ready.cns" for sample_id in config["DNA_Samples"]],
        purity="DATA/Pathological_purity_BMS_validation.txt",
        relevant_genes="DATA/TSO500_relevant_genes.txt",
        ONCOCNV_events="CNV/ONCOCNV_calls/cnv_event.txt",
        bed_file="CNV/bed/cnvkit_manifest.target.bed",
    output:
        relevant_cnvs="Results/DNA/CNV/relevant_cnv.txt",
        raw_cnv="CNV/cnvkit_calls/cnv_raw_event.txt",
    log:
        "logs/CNV_cnvkit/Filter_cnv.log",
    singularity:
        config["singularity"].get("python", config["singularity"].get("default", ""))
    script:
        "../../../scripts/python/Filter_cnv.py"
