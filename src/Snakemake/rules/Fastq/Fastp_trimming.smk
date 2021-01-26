# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf, Patrik Smeds"
__copyright__ = "Copyright 2021, Patrik Smeds, Jonas Almlöf"
__email__ = "jonas.almlöf@scilifelab.uu.se, patrik.smeds@scilifelab.uu.se"
__license__ = "GPL3"

"""
 Trims pe reads using fastq.
 Input, output and config
 ------------------------------------------------------------------------------
 Input variable: fastp_trimming_input: optional
     Default:
        "fastq_temp/{sample}_R1.fastq.gz",
        "fastq_temp/{sample}_R2.fastq.gz"
 Output variable:  
    fastp_trimming_output: optional
        Default:
            "fastq/DNA/{sample}_R1.fastq.gz",
            "fastq/DNA/{sample}_R2.fastq.gz"
    fastp_trimming_output_html: optional
        Default:
            "fastq/DNA/{sample}.html"
    fastp_trimming_output_json: optional
        Default:
            "fastq/DNA/{sample}.json"
         
 Config dict keys: values
    config["fastp"]["adapters"] optional
    config["fastp"]["extra"] optional
    config["singularity"]["fastp"] required
 Overriding input and output
 ------------------------------------------------------------------------------
 Required wildcards:
    sample
 Override input format
 Ex
  fastp_trimming_input=[
                       "fastq/{sample}.R1.fastq.gz",
                       "fastq/{sample}.R2.fastq.gz"
                       ]
 Override output format
 Ex
   fastp_trimming_output="alignment/{sample}.cutadapt.bam"
"""

_fastp_trimming_input = ["fastq_temp/DNA/{sample}_R1.fastq.gz", "fastq_temp/DNA/{sample}_R2.fastq.gz"]
try:
    _fastp_trimming_input = fastp_trimming_input
except:
    pass

_fastp_trimming_output = ["fastq/DNA/{sample}_R1.fastq.gz", "fastq/DNA/{sample}_R2.fastq.gz"]
try:
    _fastp_trimming_output = bwa_mem_output
except:
    pass

_fastp_trimming_output_html = "fastq/DNA/{sample}.html"
try:
    _fastp_trimming_output_html = fastp_trimming_output_html
except:
    pass

_fastp_trimming_output_json = "fastq/DNA/{sample}.json"
try:
    _fastp_trimming_output_json = fastp_trimming_output_json
except:
    pass


rule fastp:
    input:
        sample=_fastp_trimming_input,
    output:
        trimmed=_fastp_trimming_output,
        html=_fastp_trimming_output_html,
        json=_fastp_trimming_output_json,
    params:
        adapters=config.get("fastp", {}).get("adapters", ""),
        extra=config.get("fastp", {}).get("extra", ""),
    log:
        "logs/trimming/fastp/{sample}.log",
    threads: 5
    benchmark:
        repeat("benchmarks/trimming/fastp/{sample}.tsv", config.get("benchmark", {}).get("repeats", 1))
    singularity:
        config["singularity"].get("fastp", config["singularity"].get("default", ""))
    wrapper:
        "v0.69.0/bio/fastp"
