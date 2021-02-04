localrules: add_header_to_vcf

import src.lib.python.utils as utils

rule add_header_to_vcf:
    input:
        "{vcf}.fixChr.vcf",
    output:
        temp("{vcf}.chrAdded.vcf"),
    singularity:
        config["singularity"].get("python", config["singularity"].get("default", ""))
    params:
        type="CONTIG",
        entries=utils.create_chr_entries_for_vff_header(config['reference']['ref'] + ".fai", config['reference']['assembly'])
    script:
        "../../../scripts/python/insert_entries_vcf_header.py"