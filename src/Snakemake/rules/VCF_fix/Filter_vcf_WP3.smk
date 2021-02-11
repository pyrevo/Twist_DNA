
rule Filter_SNP_vcf:
    input:
        vcf="haplotypecaller/{sample}.VEP.vcf.gz",
    output:
        vcf="haplotypecaller/{sample}.VEP.filteredSNP.vcf.gz",
    log:
        "logs/variantCalling/Filter1/{sample}.log",
    singularity:
        config["singularity"].get("bcftools", config["singularity"].get("default", ""))
    shell:
        "bcftools filter -O v --soft-filter 'GATKCutoffSNP' "
        "-e 'TYPE=\"snp\" && "
        "(MQRankSum < -12.5 || ReadPosRankSum < -8.0 || QD < 2.0 || FS > 60.0 || "
        "(QD < 10.0 && AD[0:1] / (AD[0:1] + AD[0:0]) < 0.25 && ReadPosRankSum < 0.0) || MQ < 30.0)' "
        "-m '+' {input.vcf} "
        "| sed 's/\\\"//g' "
        "| bgzip -c > {output.vcf}"


rule Filter_INDEL_vcf:
    input:
        vcf="haplotypecaller/{sample}.VEP.filteredSNP.vcf.gz",
    output:
        vcf="haplotypecaller/{sample}.VEP.filteredSNP.filteredINDEL.vcf.gz",
    log:
        "logs/variantCalling/Filter1/{sample}.log",
    singularity:
        config["singularity"].get("bcftools", config["singularity"].get("default", ""))
    shell:
        "bcftools filter -O v --soft-filter 'GATKCutoffIndel' "
        "-e 'TYPE=\"indel\" && "
        "(ReadPosRankSum < -20.0 || QD < 2.0 || FS > 200.0 || SOR > 10.0 "
        "|| (QD < 10.0 && AD[0:1] / (AD[0:1] + AD[0:0]) < 0.25 && ReadPosRankSum < 0.0))' "
        "-m '+' {input.vcf} "
        "| sed 's/\\\"//g' "
        "| bgzip -c > {output.vcf}"
