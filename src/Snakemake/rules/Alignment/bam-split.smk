# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Patrik Smeds"
__copyright__ = "Copyright 2021, Patrik Smeds, Jonas Almlöf"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "MIT"

"""
 Rule that will split a bam file into multiple sub files, by reference.
 Input, output and config
 ------------------------------------------------------------------------------
 Input variable:
    bwa_split_input: optional
        Default:
            "alignment/{sample}.sort.bam"
 Output variable:
    bwa_split_output: optional
        Default:
            "alignment/temp/{sample}.{chr}.sort.bam"
 Config dict keys: values
    config["singularity"]["samtools"]' or config["singularity"]["default"]'  : required
 Overriding input and output
 ------------------------------------------------------------------------------
 Required wildcards:
    sample
 Override input format
 Ex
  bwa_split_input="alignment/{sample}.test.bam"
 Override output format
 Ex
   bwa_split_output="alignment/{sample}.test.{chr}.bam"
"""


_bam_split_input = "alignment/{sample}.sort.bam"
try:
    _bam_split_input = bam_split_input
except:
    pass

_bam_split_output = "alignment/temp/{sample}.{chr}.sort.bam"
try:
    _bam_split_output = bam_split_output
except:
    pass


rule bam_split:
    input:
        _bam_split_input,
    output:
        _bam_split_output
    singularity:
        config["singularity"].get("samtools", config["singularity"].get("default", ""))
    log:
        "logs/bam/split_bam_{sample}.{chr}.log",
    shell:
        "(samtools view -b {input.bam} {wildcards.chr} > {output.bam}) &> {log}"
