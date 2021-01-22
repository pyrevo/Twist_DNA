
#Todo: PoN, workflow, Twist_DNA.smk, Put approatite files into results, Filtering (only large CNV:s)?

rule collectReadCounts:
    input:
        bam="DNA_bam/{sample}-ready.bam",
        bai="DNA_bam/{sample}-ready.bam.bai",
        interval=config["CNV"]["interval"],
    output:
        "CNV_GATK/{sample}.counts.hdf5",
    params:
        mergingRule="OVERLAPPING_ONLY",
    log:
        "logs/CNV_GATK/{sample}.collectReadCounts.log"
    singularity:
        config["singularity"]["gatk4"]
    shell:
          "(gatk --java-options '-Xmx4g' CollectReadCounts -I {input.bam} -L {input.interval} \
          --interval-merging-rule {params.mergingRule} -O {output} ) &> {log}"

rule denoiseReadCounts:
    input:
        hdf5PoN=config["CNV"]["PoN"],
        hdf5Tumor="CNV_GATK/{sample}.counts.hdf5",
    output:
        denoisedCopyRatio="CNV_GATK/{sample}_clean.denoisedCR.tsv",
    params:
        stdCopyRatio="CNV_GATK/{sample}_clean.standardizedCR.tsv",
    log:
        "logs/CNV_GATK/{sample}-denoise.log"
    singularity:
        config["singularity"]["gatk4"]
    shell:
        "(gatk --java-options '-Xmx4g' DenoiseReadCounts -I {input.hdf5Tumor} \
        --count-panel-of-normals {input.hdf5PoN} \
        --standardized-copy-ratios {params.stdCopyRatio} \
        --denoised-copy-ratios {output.denoisedCopyRatio} ) &> {log}"

rule collectAllelicCounts:
    input:
        intervalList=config["CNV"]["interval"],
        bam="DNA_bam/{sample}-ready.bam",
        bai="DNA_bam/{sample}-ready.bam.bai",
        ref=config["reference"]["ref"],
    output:
        "CNV_GATK/{sample}_clean.allelicCounts.tsv",
    log:
        "logs/CNV_GATK/{sample}_allelicCounts.log"
    singularity:
        config["singularity"]["gatk4"]
    shell:
        "(gatk --java-options '-Xmx4g' CollectAllelicCounts -L {input.intervalList} \
        -I {input.bam} -R {input.ref} \
        -O {output} ) &> {log}"


rule modelSegments:
    input:
        denoisedCopyRatio="CNV_GATK/{sample}_clean.denoisedCR.tsv",
        allelicCounts="CNV_GATK/{sample}_clean.allelicCounts.tsv"
    output:
        #"CNV_GATK/{sample}_clean.modelBegin.seg",
        "CNV_GATK/{sample}_clean.modelFinal.seg",
        "CNV_GATK/{sample}_clean.cr.seg",
        #"CNV_GATK/{sample}_clean.modelBegin.af.param",
        #"CNV_GATK/{sample}_clean.modelBegin.cr.param",
        #"CNV_GATK/{sample}_clean.modelFinal.af.param",
        #"CNV_GATK/{sample}_clean.modelFinal.cr.param",
        "CNV_GATK/{sample}_clean.hets.tsv",
    params:
        outDir="CNV_GATK/",
        outPrefix="{sample}_clean",
    log:
        "logs/CNV_GATK/{sample}_modelSegments.log"
    singularity:
        config["singularity"]["gatk4"]
    shell:
        "(gatk --java-options '-Xmx4g' ModelSegments \
        --denoised-copy-ratios {input.denoisedCopyRatio} \
        --allelic-counts {input.allelicCounts} \
        --output {params.outDir} --output-prefix {params.outPrefix} ) &> {log}"

rule callCopyRatioSegments:
    input:
        "CNV_GATK/{sample}_clean.cr.seg",
    output:
        "CNV_GATK/{sample}_clean.calledCNVs.seg",
    log:
        "logs/CNV_GATK/{sample}_calledCRSegments.log"
    singularity:
        config["singularity"]["gatk4"]
    shell:
        "(gatk CallCopyRatioSegments --input {input} \
        --output {output} ) &> {log}"

rule plotModeledSegments:
    input:
        denoisedCopyRatio="CNV_GATK/{sample}_clean.denoisedCR.tsv",
        allelicCounts="CNV_GATK/{sample}_clean.hets.tsv",
        segments="CNV_GATK/{sample}_clean.modelFinal.seg",
        refDict=config["reference"]["ref"][:-5]+"dict",
    output:
        "CNV_GATK/{sample}_clean.calledCNVs.modeled.png",
    params:
        outDir="CNV_GATK/",
        outPrefix="{sample}_clean.calledCNVs",
        pointSize=2.0,
    log:
        "logs/CNV_GATK/{sample}_plotSegments.log"
    singularity:
        config["singularity"]["gatk4"]
    shell:
        "(gatk PlotModeledSegments --denoised-copy-ratios {input.denoisedCopyRatio} \
        --allelic-counts {input.allelicCounts} --segments {input.segments} \
        --sequence-dictionary {input.refDict} \
        --point-size-allele-fraction {params.pointSize} --point-size-copy-ratio {params.pointSize} \
        --output {params.outDir} --output-prefix {params.outPrefix} ) &> {log} "
