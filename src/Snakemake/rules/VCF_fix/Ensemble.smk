
rule ensemble_filter:
    input:
        vcf="recall/{sample}.ensemble.vcf.gz",
    output:
        vcf="recall/{sample}.ensemble.final.vcf.gz",
    shell:
        "python3 src/filter_by_num_callers.py -v {input.vcf} -d | bgzip > {output.vcf} && "
        "tabix {output.vcf}"
