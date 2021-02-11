
rule Extract_chrX_vcf:
    input:
        vcf="haplotypecaller/{sample}.vep.filteredSNP.filteredINDEL.filteredAF.Cartagenia.noHLA.vcf",
        bed=config["bed"]["bedfile_chrX"]
    output:
        vcf="haplotypecaller/{sample}.vep.filteredSNP.filteredINDEL.filteredAF.Cartagenia.noHLA.chrX.vcf",
    singularity:
        config["singularity"].get("bedtools", config["singularity"].get("default", ""))
    shell:
        "bedtools intersect -header -a {input.vcf} -b {input.bed} > {output.vcf}"


rule chrX_vcfstats:
    input:
        vcf="haplotypecaller/{sample}.vep.filteredSNP.filteredINDEL.filteredAF.Cartagenia.noHLA.chrX.vcf",
    output:
        vcf="haplotypecaller/{sample}.vep.filteredSNP.filteredINDEL.filteredAF.Cartagenia.noHLA.chrX.vcfstats.vcf",
    singularity:
        config["singularity"].get("vcftools", config["singularity"].get("default", ""))
    shell:
        "vcf-stats {input.vcf} > {output.vcf}"
