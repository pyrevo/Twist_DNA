
rule Haplotypecaller:
    input:
        bam="alignment/{sample}.bam",
    output:
        vcf="haplotypecaller/{sample}.haplotypecaller.fixAF.vcf.gz",
    params:
        reference=config["reference"]["ref"],
        bed=config["bed"]["bedfile"],
        annotation="--annotation MappingQualityRankSumTest --annotation MappingQualityZero --annotation QualByDepth \
                    --annotation ReadPosRankSumTest --annotation RMSMappingQuality --annotation BaseQualityRankSumTest \
                    --annotation FisherStrand --annotation MappingQuality --annotation DepthPerAlleleBySample \
                    --annotation Coverage --annotation ClippingRankSumTest --annotation DepthPerSampleHC",
        extra="--interval-set-rule INTERSECTION --native-pair-hmm-threads 1 -ploidy 2",
    log:
    singularity:
        config["singularity"].get("gatk4", config["singularity"].get("default", ""))
    shell:
        "gatk --java-options '-Xmx6g' HaplotypeCaller -R {params.reference} {params.annotation} -I {input.bam} "
        "-L {params.bed} {params.extra} --output {output.gz}"
