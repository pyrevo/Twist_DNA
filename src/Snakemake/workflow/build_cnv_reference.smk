
include: "../rules/Alignment/index_bam.smk"
include: "../rules/Fastq/fix_fastq_DNA.smk"

if config["programs"]["Trimming"] == "Cutadapt":

    include: "../rules/Fastq/Cutadapt_trimming.smk"

elif config["programs"]["Trimming"] == "Fastp":

    include: "../rules/Fastq/Fastp_trimming.smk"

else:

    include: "../rules/Fastq/move_fastq.smk"

if config["programs"]["Duplicates"] == "fgbio":

    bam_split_input = "Bam/DNA/{sample}-ready.bam"

    include: "../rules/Alignment/fgbio.smk"


else:
    if config["programs"]["markduplicate"] == "GPU":

        include: "../rules/Alignment/GPU_alignment.smk"


    elif config["programs"]["markduplicate"] == "picard":

        bam_split_input = "alignment/{sample}.sort.bam"
        markduplicates_input = "alignment/{sample}.{chr}.bam"
        markduplicates_output = "alignment/{sample}.dedup.{chr}.bam"
        bam_merge_input = "alignment/{sample}.dedup.__CHR__.bam"
        bam_merge_output = "Bam/DNA/{sample}-ready.bam"

        include: "../rules/Alignment/bwa-mem.smk"
        include: "../rules/Alignment/MarkDuplicates.smk"


    elif config["programs"]["markduplicate"] == "picard_UMI":

        bwa_mem_output = "alignment/{sample}.sort.noUMI.bam"
        #bam_split_input = "alignment/{sample}.sort.bam"
        #markduplicates_input = "alignment/{sample}.{chr}.bam"
        markduplicates_output = "alignment/{sample}.dedup.bam"
        #bam_merge_input = "alignment/{sample}.dedup.__CHR__.bam"
        #bam_merge_output = "Bam/DNA/{sample}-ready.bam"
        #mutect_input = "alignment/{sample}.dedup.{chr}.bam"
        #freebayes_input = "alignment/{sample}.dedup.{chr}.bam"
        #vardict_input = "alignment/{sample}.dedup.{chr}.bam"

        include: "../rules/Alignment/bwa-mem.smk"
        include: "../rules/Alignment/MarkDuplicatesUMI.smk"


include: "../rules/Alignment/bam-split.smk"
include: "../rules/Alignment/bam-merge.smk"

include: "../rules/CNV/cnvkit_PoN.smk"
include: "../rules/CNV/GATK_CNV_PoN.smk"
