

rule intron_filter:
    input:
        vcf="recall/{sample}.ensemble.vep.vcf.gz",
        bed=config["bed"]["bedfile"],
    output:
        vcf="recall/{sample}.ensemble.vep.exon.vcf.gz",
    params:
        vcf="recall/{sample}.ensemble.vep.exon.vcf",
    # singularity:
    #    config["singularity"]["python"]
    singularity:
        config["singularity"].get("python_htslib", config["singularity"].get("default", ""))
    shell:
        "python src/scripts/python/filter_introns.py {input.vcf} {input.bed} {params.vcf} &&"
        "bgzip {params.vcf} && "
        "tabix {output.vcf}"


rule Soft_filter:
    input:
        vcf="recall/{sample}.ensemble.vep.exon.vcf.gz",
    output:
        vcf="recall/{sample}.ensemble.vep.exon.soft_filter.vcf",
    params:
        filter="-e 'FORMAT/AD<20 || FORMAT/DP<50 || FORMAT/AF<0.05'",
    singularity:
        config["singularity"].get("bcftools", config["singularity"].get("default", ""))
        # "/projects/wp2/nobackup/Twist_Myeloid/Containers/bcftools-1.9--8.simg"
    shell:
        "bcftools filter -O v -o {output.vcf} --soft-filter 'Soft_filter' {params.filter} -m '+' {input.vcf}"


rule ffpe_filter:
    input:
        vcf="recall/{sample}.ensemble.vep.exon.soft_filter.vcf",
        bam="DNA_bam/{sample}-ready.bam",
        bai="DNA_bam/{sample}-ready.bam.bai",
    params:
        vcf_ffpe_temp=temp("recall/{sample}.ensemble.vep.exon.soft_filter.ffpe.temp.vcf"),
        vcf_ffpe=temp("recall/{sample}.ensemble.vep.exon.soft_filter.ffpe.vcf"),
        java=config["java"]["SOBDetector"],
    output:
        gvcf="recall/{sample}.ensemble.vep.exon.soft_filter.vcf.gz",
        gvcf_ffpe="recall/{sample}.ensemble.vep.exon.soft_filter.ffpe.vcf.gz",
    shell:
        #"module load oracle-jdk-1.8/1.8.0_162 && "
        "java -jar {params.java} --input-type VCF --input-variants {input.vcf} --input-bam {input.bam} --output-variants {params.vcf_ffpe_temp} && "
        "python src/scripts/python/Add_FFPE_column_to_vcf.py {params.vcf_ffpe_temp} {params.vcf_ffpe} && "
        "bgzip {params.vcf_ffpe} && "
        "tabix {output.gvcf_ffpe} && "
        "bgzip {input.vcf} && "
        "tabix {output.gvcf}"
