

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
    shell:
        "bcftools filter -O v -o {output.vcf} --soft-filter 'Soft_filter' {params.filter} -m '+' {input.vcf}"


rule ffpe_filter:
    input:
        vcf="recall/{sample}.ensemble.vep.exon.soft_filter.vcf",
        bam="Bam/DNA/{sample}-ready.bam",
        bai="Bam/DNA/{sample}-ready.bam.bai",
    params:
        vcf_ffpe_temp=temp("Results/DNA/{sample}/vcf/{sample}.ensemble.vep.exon.soft_filter.ffpe.temp.vcf"),
        vcf_ffpe=temp("Results/DNA/{sample}/vcf/{sample}.ensemble.vep.exon.soft_filter.ffpe.vcf"),
        java=config["java"]["SOBDetector"],
    output:
        vcf_gz="Results/DNA/{sample}/vcf/{sample}.ensemble.vep.exon.soft_filter.vcf.gz",
        vcf_gz_ffpe="Results/DNA/{sample}/vcf/{sample}.ensemble.vep.exon.soft_filter.ffpe.vcf.gz",
    shell:
        #"module load oracle-jdk-1.8/1.8.0_162 && "
        "java -jar {params.java} --input-type VCF --input-variants {input.vcf} --input-bam {input.bam} --output-variants {params.vcf_ffpe_temp} && "
        "python src/scripts/python/Add_FFPE_column_to_vcf.py {params.vcf_ffpe_temp} {params.vcf_ffpe} && "
        "bgzip {params.vcf_ffpe} && "
        "tabix {output.vcf_gz_ffpe} && "
        "bgzip -c {input.vcf} > {output.vcf_gz} && "
        "tabix {output.vcf_gz}"
