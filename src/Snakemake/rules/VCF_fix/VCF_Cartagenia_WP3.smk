
# rule gunzip_VCF_Cartagenia:
#     input:
#         vcf="haplotypecaller/{sample}.vep.filteredSNP.filteredINDEL.filteredAF.vcf.gz",
#     output:
#         vcf="haplotypecaller/{sample}.vep.filteredSNP.filteredINDEL.filteredAF.vcf",
#     shell:
#         "gunzip -c {input.vcf} > {output.vcf}"


rule VCF_Cartagenia:
    input:
        vcf="haplotypecaller/{sample}.vep.filteredSNP.filteredINDEL.filteredAF.vcf",
    output:
        vcf="haplotypecaller/{sample}.vep.filteredSNP.filteredINDEL.filteredAF.Cartagenia.vcf",
    run:
        import subprocess
        subprocess.call("src/scripts/perl/Vcf_to_Cartagenia.pl " + input.vcf + " " + output.vcf, shell=True)


rule gzip_VCF_Cartagenia:
    input:
        vcf="haplotypecaller/{sample}.vep.filteredSNP.filteredINDEL.filteredAF.Cartagenia.vcf",
    output:
        vcf="haplotypecaller/{sample}.vep.filteredSNP.filteredINDEL.filteredAF.Cartagenia.vcf.gz",
    shell:
        "gzip {input.vcf}"
