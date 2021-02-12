#
# rule Haplotypecaller:
#     input:
#         bam="alignment/{sample}.bam",
#         bai="alignment/{sample}.bam.bai",
#     output:
#         vcf="haplotypecaller/{sample}.vcf.gz",
#     params:
#         reference=config["reference"]["ref"],
#         bed=config["bed"]["bedfile"],
#         annotation="--annotation MappingQualityRankSumTest --annotation MappingQualityZero --annotation QualByDepth \
#                     --annotation ReadPosRankSumTest --annotation RMSMappingQuality --annotation BaseQualityRankSumTest \
#                     --annotation FisherStrand --annotation MappingQuality --annotation DepthPerAlleleBySample \
#                     --annotation Coverage --annotation ClippingRankSumTest --annotation DepthPerSampleHC",
#         #extra="--interval-set-rule INTERSECTION --native-pair-hmm-threads 1 -ploidy 2",
#         #Add padding?
#         extra="--native-pair-hmm-threads 1 -ploidy 2",
#     log:
#         "logs/variantCalling/Haplotypecaller/call/{sample}.log",
#     singularity:
#         config["singularity"].get("gatk4", config["singularity"].get("default", ""))
#     shell:
#         "(gatk --java-options '-Xmx6g' HaplotypeCaller -R {params.reference} {params.annotation} -I {input.bam} "
#         "-L {params.bed} {params.extra} --output {output.vcf}) &> {log}"


chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13",
               "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]


rule Split_bed:
    input:
        bed=config["bed"]["bedfile"],
    output:
        beds=expand("bedfiles/Twist_Exome_Target_chr{chrom}.bed", chrom=chromosomes),
    log:
        "logs/variantCalling/Haplotypecaller/call/{sample}.{chrom}.log",
    run:
        import subprocess
        for chrom in chromosomes :
            subprocess.call("grep -P \"^" + chrom + "\t\" " + input.bed + " > bedfiles/Twist_Exome_Target_chr" + chrom + ".bed", shell=True)

rule Haplotypecaller:
    input:
        bam="alignment/{sample}.bam",
        bai="alignment/{sample}.bam.bai",
        bed="bedfiles/Twist_Exome_Target_chr{chrom}.bed"
    output:
        vcf=temp("haplotypecaller/{sample}.{chrom}.vcf.gz"),
    params:
        reference=config["reference"]["ref"],
        annotation="--annotation MappingQualityRankSumTest --annotation MappingQualityZero --annotation QualByDepth \
                    --annotation ReadPosRankSumTest --annotation RMSMappingQuality --annotation BaseQualityRankSumTest \
                    --annotation FisherStrand --annotation MappingQuality --annotation DepthPerAlleleBySample \
                    --annotation Coverage --annotation ClippingRankSumTest --annotation DepthPerSampleHC",
        #extra="--interval-set-rule INTERSECTION --native-pair-hmm-threads 1 -ploidy 2",
        #Add padding?
        extra="--interval-set-rule INTERSECTION --native-pair-hmm-threads 1 -ploidy 2",
    log:
        "logs/variantCalling/Haplotypecaller/call/{sample}.{chrom}.log",
    singularity:
        config["singularity"].get("gatk4", config["singularity"].get("default", ""))
    shell:
        "(gatk --java-options '-Xmx6g' HaplotypeCaller -R {params.reference} {params.annotation} -I {input.bam} "
        "-L {input.bed} {params.extra} --output {output.vcf}) &> {log}"


rule Merge_Haplotypecaller_vcf:
    input:
        vcf=expand("haplotypecaller/{{sample}}.{chrom}.vcf.gz", chrom=chromosomes),
    output:
        vcf="haplotypecaller/{sample}.vcf.gz"
    log:
        "logs/variantCalling/Haplotypecaller/merge_vcf/{sample}.log",
    singularity:
        config["singularity"].get("bcftools", config["singularity"].get("default", ""))
    wrapper:
        "0.70.0/bio/bcftools/concat"
