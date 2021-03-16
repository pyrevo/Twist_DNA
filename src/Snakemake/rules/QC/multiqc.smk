
rule multiqcBatch:
    input:
        qc1=expand("qc/{sample}/{sample}_R1_fastqc.zip", sample=config["DNA_Samples"]),
        qc2=expand("qc/{sample}/{sample}_R2_fastqc.zip", sample=config["DNA_Samples"]),
        qc3=expand("qc/{sample}/{sample}.samtools-stats.txt", sample=config["DNA_Samples"]),
        qc4=expand("qc/{sample}/{sample}.HsMetrics.txt", sample=config["DNA_Samples"]),
        qc5=expand("qc/{sample}/{sample}_batchStats.done", sample=config["DNA_Samples"]),  #Wait until all in table
        qc6="qc/batchQC_stats_mqc.json",
        #qc6 = expand("qc/{sample}/{sample}_avg_CV_genes_over_500X.txt", sample=config["DNA_Samples"]),
        #qc7=expand("qc/{sample}/{sample}.gc_bias.summary_metrics.txt", sample=config["DNA_Samples"]),
        qc7=expand("qc/{sample}/{sample}.duplication_metrics.txt", sample=config["DNA_Samples"]),
    output:
        "qc/MultiQC.html",
    params:
        "-c " + config["configfiles"]["multiqc"] + " --ignore *_stats_mqc.csv",
    log:
        "logs/report/multiqc.log",
    singularity:
        config["singularity"].get("multiqc", config["singularity"].get("default", ""))
    wrapper:
        "0.72.0/bio/multiqc"


rule copy_multiqc:
    input:
        "qc/MultiQC.html",
    output:
        "Results/DNA/MultiQC.html",
    log:
        "logs/report/multiqc_copy.log",
    shell:
        "(cp qc/MultiQC.html Results/DNA/MultiQC.html) &> {log}"
