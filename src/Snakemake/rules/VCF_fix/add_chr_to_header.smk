# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Patrik Smeds"
__copyright__ = "Copyright 2021, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "GPL3"

"""
 
 Rule used to add a contigs to a vcf file. Uses fai file to create contig information.
 
 """


localrules:
    add_header_to_vcf,


import src.lib.python.utils as utils


rule add_header_to_vcf:
    input:
        "{vcf}.fixChr.vcf",
    output:
        temp("{vcf}.chrAdded.vcf"),
    singularity:
        config["singularity"].get("python", config["singularity"].get("default", ""))
    params:
        type="contig",
        entries=utils.create_chr_entries_for_vff_header(
            config['reference']['ref'] + ".fai", config['reference']['assembly']
        ),
    script:
        "../../../scripts/python/insert_entries_vcf_header.py"
