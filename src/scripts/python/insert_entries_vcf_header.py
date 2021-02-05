#!/bin/python3
from pysam import VariantFile

vcf_in = VariantFile(snakemake.input[0])
new_header = vcf_in.header

for entry in snakemake.params.entries:
    new_header.add_meta(snakemake.params.type, items=entry)

vcf_out = VariantFile(snakemake.output[0], 'w', header=new_header)
for record in vcf_in.fetch():
    vcf_out.write(record)
