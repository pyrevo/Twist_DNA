
rule VCF_Cartagenia:
    input:
        vcf="haplotypecaller/{sample}.vep.filteredSNP.filteredINDEL.vcf.gz",
    output:
        vcf="haplotypecaller/{sample}.vep.filteredSNP.filteredINDEL.Cartagenia.vcf.gz",
    run:
        import subprocess
        subprocess.call("src/scripts/perl/Vcf_bcbio_to_Cartagenia.pl " + input.vcf + " " + output.vcf, shell=True)
