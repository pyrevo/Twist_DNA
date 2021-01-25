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
            "alignment/{sample}.bam"
 Output variable:   
    bwa_split_output: optional
        Default: where REF_ will be converted to REF_chr1, ... REF_chrX if chromosome 1 is named chr1 in the fai file
            ["alignment/{sample}.REF_.bam", ... , "alignment/{sample}.REF_.bam"]    
 Config dict keys: values
    config["reference"]["ref"]': required for default output, used to locate fai-file
    config["singularity"]["bamtools"]' or config["singularity"]["default"]'  : required 
 Overriding input and output
 ------------------------------------------------------------------------------
 Required wildcards:
    sample
 Override input format
 Ex
  bwa_split_input="alignment/{sample}.test.bam"
 Override output format
 Ex
   bwa_split_output="alignment/{sample}.test.REF_.bam"
"""


def extract_chr(file):
    chr = None
    with open(file) as lines:
        chr = [line.split("\t")[0] for line in lines]
    print(str(chr))
    return chr


_bam_split_input = "alignment/{sample}.bam"
try:
    _bam_split_input = bam_split_input
except:
    pass

_bam_split_output = "alignment/{sample}.REF_.bam"
try:
    _bam_split_output = bam_split_output
except:
    pass


rule bam_split:
    input:
        _bam_split_input,
    output:
        [_bam_split_output.replace("REF_", "REF_" + chr) for chr in extract_chr(config['reference']['ref'] + ".fai")],
    singularity:
        config["singularity"].get("bamtools", config["singularity"].get("default", ""))
    wrapper:
        "file:///projects/wp4/nobackup/workspace/snakemake-wrappers/bio/bamtools/split"
