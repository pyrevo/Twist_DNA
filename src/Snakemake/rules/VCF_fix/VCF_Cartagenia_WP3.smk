
rule gunzip_VCF_Cartagenia:
    input:
        vcf="haplotypecaller/{sample}.vep.filteredSNP.filteredINDEL.vcf.gz",
    output:
        vcf="haplotypecaller/{sample}.vep.filteredSNP.filteredINDEL.vcf",
    shell:
        "gunzip -c {input.vcf} > {output.vcf}"


rule VCF_Cartagenia:
    input:
        vcf="haplotypecaller/{sample}.vep.filteredSNP.filteredINDEL.vcf",
    output:
        vcf="haplotypecaller/{sample}.vep.filteredSNP.filteredINDEL.Cartagenia.vcf",
    run:
        import subprocess
        subprocess.call("src/scripts/perl/Vcf_to_Cartagenia.pl " + input.vcf + " " + output.vcf, shell=True)


rule gzip_VCF_Cartagenia:
    input:
        vcf="haplotypecaller/{sample}.vep.filteredSNP.filteredINDEL.Cartagenia.vcf",
    output:
        vcf="haplotypecaller/{sample}.vep.filteredSNP.filteredINDEL.Cartagenia.vcf.gz",
    shell:
        "gzip {input.vcf}"
