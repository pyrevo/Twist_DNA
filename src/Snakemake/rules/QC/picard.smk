localrules:
    touchBatch,


rule picardHsMetrics:
    input:
        bam="DNA_bam/{sample}-ready.bam",
        intervals=config["bed"]["intervals"],
    output:
        "qc/{sample}/{sample}.HsMetrics.txt",
    log:
        "logs/qc/picard/HsMetrics/{sample}.log",
    singularity:
        config["singularity"]["picard"]
    shell:
        "(java -Xmx4g -jar /opt/conda/share/picard-2.20.1-0/picard.jar CollectHsMetrics BAIT_INTERVALS={input.intervals} TARGET_INTERVALS={input.intervals} INPUT={input.bam} OUTPUT={output}) &> {log}"


rule picardInsertSize:
    input:
        bam="DNA_bam/{sample}-ready.bam",
    output:
        txt="qc/{sample}/{sample}.insert_size_metrics.txt",
        pdf="qc/{sample}/{sample}.insert_size_histogram.pdf",
    log:
        "logs/qc/picard/InsertSize/{sample}.log",
    singularity:
        config["singularity"]["picard"]
    shell:
        "(java -Xmx4g -jar /opt/conda/share/picard-2.20.1-0/picard.jar CollectInsertSizeMetrics INPUT={input.bam} O={output.txt} H={output.pdf}) &> {log}"


rule PicardAlignmentSummaryMetrics:
    input:
        bam="DNA_bam/{sample}-ready.bam",
        ref=config["reference"]["ref"],
    output:
        "qc/{sample}/{sample}.alignment_summary_metrics.txt",
    log:
        "logs/qc/picard/AlignmentSummaryMetrics/{sample}.log",
    singularity:
        config["singularity"]["picard"]
    shell:
        "(java -Xmx4g -jar /opt/conda/share/picard-2.20.1-0/picard.jar CollectAlignmentSummaryMetrics INPUT={input.bam} R={input.ref} OUTPUT={output}) &> {log}"


rule GcBiasSummaryMetrics:
    input:
        bam="DNA_bam/{sample}-ready.bam",
        ref=config["reference"]["ref"],
    output:
        summary="qc/{sample}/{sample}.gc_bias.summary_metrics.txt",
        gc="qc/{sample}/{sample}.gc_bias.metrics.txt",
        pdf="qc/{sample}/{sample}.gc_bias.metrics.pdf",
    log:
        "logs/qc/picard/GcBiasSummaryMetrics/{sample}.log",
    singularity:
        config["singularity"]["picard"]
    shell:
        "(java -Xmx4g -jar /opt/conda/share/picard-2.20.1-0/picard.jar CollectGcBiasMetrics I={input.bam} R={input.ref} O={output.gc} CHART={output.pdf} S={output.summary}) &> {log}"
