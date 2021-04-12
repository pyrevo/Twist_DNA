# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Patrik Smeds, Jonas Almlöf"
__copyright__ = "Copyright 2021, Patrik Smeds, Jonas Almlöf"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "MIT"

"""
 Rule that will mark duplicates using UMIs

 """


_markduplicates_input = "alignment/{sample}.{chr}.bam"
try:
    _markduplicates_input = markduplicates_input
except:
    pass


_markduplicates_output = "alignment/{sample}.dup.{chr}.bam"
try:
    _markduplicates_output = markduplicates_output
except:
    pass


# rule UmiAwareMarkDuplicatesWithMateCigar:
#     input:
#         bam=_markduplicates_input,
#         bai=_markduplicates_input + ".bai",
#     output:
#         bam=temp(_markduplicates_output),
#     params:
#         metric=temp("alignment/{sample}_DuplicationMetrics.{chr}.txt"),
#         UMI_metric="alignment/{sample}_UmiMetrics.{chr}.txt",
#     log:
#         "logs/map/MarkDup/{sample}-ready.{chr}.log",
#     threads: 2
#     container:
#         config["singularity"].get("picard", config["singularity"].get("default", ""))
#     shell:
#         "(picard UmiAwareMarkDuplicatesWithMateCigar DUPLEX_UMI=true MAX_EDIT_DISTANCE_TO_JOIN=1 BARCODE_TAG=RX INPUT={input.bam} OUTPUT={output.bam} METRICS_FILE={params.metric} UMI_METRICS={params.UMI_metric}) &> {log}"
#

rule UmiMarkDuplicates:
    input:
        bam=_markduplicates_input,
        bai=_markduplicates_input + ".bai",
    output:
        bam=temp(_markduplicates_output),
    params:
        metric=temp("alignment/{sample}_DuplicationMetrics.{chr}.txt"),
    log:
        "logs/map/MarkDup/{sample}-ready.{chr}.log",
    threads: 2
    container:
        config["singularity"].get("picard", config["singularity"].get("default", ""))
    shell:
        "(picard MarkDuplicates DUPLEX_UMI=true BARCODE_TAG=RX INPUT={input.bam} OUTPUT={output.bam} METRICS_FILE={params.metric}) &> {log}"
