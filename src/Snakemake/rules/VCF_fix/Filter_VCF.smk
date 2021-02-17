

rule intron_filter:
    input:
        vcf="recall/{sample}.ensemble.final.vcf.gz",
        bed=config["bed"]["bedfile"],
    output:
        vcf="recall/{sample}.ensemble.final.exon.vcf.gz",
    params:
        vcf="recall/{sample}.ensemble.final.exon.vcf",
    shell:
        "python3 src/filter_introns.py {input.vcf} {input.bed} {params.vcf} &&"
        "bgzip {params.vcf} && "
        "tabix {output.vcf}"


rule Soft_filter:
    input:
        vcf="recall/{sample}.ensemble.final.exon.vcf.gz",
    output:
        vcf="recall/{sample}.ensemble.final.exon.soft_filter.vcf",
    params:
        filter="-e 'FORMAT/AD<20 || FORMAT/DP<50 || FORMAT/AF<0.05'",
    singularity:
        "/projects/wp2/nobackup/Twist_Myeloid/Containers/bcftools-1.9--8.simg"
    shell:
        "bcftools filter -O v -o {output.vcf} --soft-filter 'Soft_filter' {params.filter} -m '+' {input.vcf}"


rule ffpe_filter:
    input:
        vcf="recall/{sample}.ensemble.final.exon.soft_filter.vcf",
        bam="DNA_bam/{sample}-ready.bam",
        bai="DNA_bam/{sample}-ready.bam.bai",
    params:
        vcf_ffpe_temp=temp("recall/{sample}.ensemble.final.exon.soft_filter.ffpe.temp.vcf"),
        vcf_ffpe=temp("recall/{sample}.ensemble.final.exon.soft_filter.ffpe.vcf"),
    output:
        gvcf="recall/{sample}.ensemble.final.exon.soft_filter.vcf.gz",
        gvcf_ffpe="recall/{sample}.ensemble.final.exon.soft_filter.ffpe.vcf.gz",
    shell:
        #"module load oracle-jdk-1.8/1.8.0_162 && "
        "java -jar SOBDetector/SOBDetector_v1.0.1.jar --input-type VCF --input-variants {input.vcf} --input-bam {input.bam} --output-variants {params.vcf_ffpe_temp} && "
        "python src/Add_FFPE_column_to_vcf.py {params.vcf_ffpe_temp} {params.vcf_ffpe} && "
        "bgzip {params.vcf_ffpe} && "
        "tabix {output.gvcf_ffpe} && "
        "bgzip {input.vcf} && "
        "tabix {output.gvcf}"
