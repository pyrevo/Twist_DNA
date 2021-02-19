#!/bin/python3.6
import sys
from pysam import VariantFile
import subprocess


vcf_in = VariantFile(snakemake.input[0])
new_header = new_header = vcf_in.header
vcf_out = VariantFile(snakemake.output[0], 'w', header=new_header)

for record in vcf_in.fetch():
    if len(record.ref) != len(record.alts[0]):  # if InDel
        if (
            "mutect2" in record.info["CALLERS"] or "vardict" in record.info["CALLERS"]
        ):
            vcf_out.write(record)
    elif len(record.info["CALLERS"]) >= 2:
        vcf_out.write(record)
