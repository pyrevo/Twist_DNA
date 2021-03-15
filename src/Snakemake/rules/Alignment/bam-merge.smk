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
    bwa_merge_input: optional
        Default:
            "alignment/temp/{sample}.__CHR__.sort.bam"
        where __REF__ will be converted to chr1, ... chrX if chromosome 1 is named chr1 in the fai file
            ["alignment/{sample}.chr1.bam", ... , "alignment/{sample}.chrXX.bam"]. __CHR__ can be change by
            setting "_merge_separator"
 Output variable:
    bwa_split_output: optional
        Default:
            "alignment/{sample}.sort.bam"
 Config dict keys: values
    config["reference"]["ref"]': required for default output, used to locate fai-file
    config["singularity"]["samtools"]' or config["singularity"]["default"]'  : required
 Overriding input and output
 ------------------------------------------------------------------------------
 Required wildcards:
    sample
 Override input format
 Ex
  bwa_merge_input="alignment/{sample}.test.__CHR__.bam"
 Override output format
 Ex
   bwa_merge_output="alignment/{sample}.test.bam"
"""

import src.lib.python.utils as utils


_merge_separator = "__CHR__"
try:
    _merge_separator = merge_separator
except:
    pass


_bam_merge_input = "alignment/temp/{sample}." + _merge_separator + ".sort.bam"
try:
    _bam_merge_input = bam_merge_input
except:
    pass

_bam_merge_output = "alignment/{sample}.sort.bam"
try:
    _bam_merge_output = bam_merge_output
except:
    pass


rule bam_split:
    input:
        [_bam_split_input.replace(_merge_separator, chr)
        for chr in utils.extract_chr(config['reference']['ref'] + ".fai", filter_out=config.get("skip_chrs",['chrM'])],
    output:
        _bam_split_output
    singularity:
        config["singularity"].get("samtools", config["singularity"].get("default", ""))
    log:
        "logs/bam/merge_bam_{sample}.{chr}.log",
    shell:
        "(samtools view -b {input.bam} {wildcards.chr} > {output.bam}) &> {log}"
