# snakemake -p -j 64 --drmaa "-A wp1 -p core -n {cluster.n} -t {cluster.time}"  -s ./src/Snakemake/rules/CNV/GATK_CNV_PoN.smk --use-singularity --singularity-args "--bind /data --bind /projects --bind /scratch " --cluster-config Config/Slurm/cluster.json


configfile: "Twist_DNA_PoN.yaml"


rule all:
    input:
        "Normals/GATK4/readCountPoN.hdf5",
        expand("Normals/GATK4/{normal}.counts.hdf5", normal=config["DNA_Samples"]),


## Create interval list and preprocess interval list, only needed when updating bedfile
# rule bedToIntervalList:
#     input:
#         bed=config["bed"]["bedfile"],  ## Annotated clostest bedfile? Really needed?
#         refDict=config["reference"]["ref"],  ##Have to be a .dict in same folder as .fasta
#     output:
#         "bedFiles/TM_TE-annotated_closest-noduplicates.interval_list",  #Should be based on bedfile...
#     log:
#         "logs/Normals/TM_TE-annotated_closest-noduplicates.log",
#     singularity:
#         config["singularity"]["gatk4"]
#     shell:
#         "(gatk BedToIntervalList  -I {input.bed} -O {output} -SD {input.refDict} ) &> {log} "


rule preprocessIntervals:
    input:
        ref=config["reference"]["ref"],
        intervalList=config["bed"]["intervals"],  #targets_C.interval_list interval list picard style
    output:
        "bedFiles/pool1_pool2_nochr_3c.annotated.preprocessed.interval_list",
    params:
        binLength=0,  #WGS 1000
        mergingRule="OVERLAPPING_ONLY",
    log:
        "logs/Normals/GATK/preprocessIntervals.log",
    singularity:
        config["singularity"]["gatk4"]
    shell:
        "(gatk --java-options '-Xmx4g' PreprocessIntervals -L {input.intervalList} -R {input.ref} "
        "--bin-length {params.binLength} --interval-merging-rule {params.mergingRule} "
        "-O {output}) &> {log}"


# From here need to be redone when added new samples
rule collectReadCounts:
    input:
        #bam=lambda wildcards: config["normal"][wildcards.normal],
        bams="DNA_bam/{normal}-ready.bam",
        interval="bedFiles/pool1_pool2_nochr_3c.annotated.preprocessed.interval_list",
    output:
        "Normals/GATK4/{normal}.counts.hdf5",  #Should have date in it?
    params:
        mergingRule="OVERLAPPING_ONLY",
    log:
        "logs/Normals/GATK/{normal}.collectReadCounts.log",
    singularity:
        config["singularity"]["gatk4"]
    shell:
        "(gatk --java-options '-Xmx4g' CollectReadCounts -I {input.bam} -L {input.interval} "
        "--interval-merging-rule {params.mergingRule} -O {output}) &> {log}"


rule createReadCountPanelOfNormals:
    input:
        expand("Normals/GATK4/{normal}.counts.hdf5", normal=config["normal"]),
    output:
        "Normals/GATK4/readCountPoN.hdf5",
    params:
        minIntervalMedianPerc=5.0,
        input=lambda wildcards, input: " -I ".join(input),
    log:
        "logs/Normals/GATK/readCountPoN.log",
    singularity:
        config["singularity"]["gatk4"]
    shell:
        "(gatk --java-options '-Xmx4g' CreateReadCountPanelOfNormals -I {params.input} "
        "--minimum-interval-median-percentile {params.minIntervalMedianPerc} "
        "-O {output}) &> {log}"
