
_markduplicates_input = "alignment/temp/{sample}.sort.{chr}.bam"
try:
    _markduplicates_input = markduplicates_input
except:
    pass


_markduplicates_output = "alignment/temp/{sample}.dup.{chr}.bam"
try:
    _markduplicates_output = markduplicates_output
except:
    pass


rule MarkDuplicates:
    input:
        bam=_markduplicates_input,
        bai=_markduplicates_input + ".bai",
    output:
        bam=temp(_markduplicates_output),
    params:
        metric="qc/{sample}_DuplicationMetrics.{chr}.txt",
    log:
        "logs/map/MarkDup/{sample}-ready.{chr}.log",
    threads: 2
    container:
        config["singularity"].get("picard", config["singularity"].get("default", ""))
    shell:
        "(picard MarkDuplicates INPUT={input.bam} OUTPUT={output.bam} METRICS_FILE={params.metric}) &> {log}"
