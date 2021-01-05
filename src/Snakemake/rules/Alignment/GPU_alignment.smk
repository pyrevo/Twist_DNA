

rule GPU_align:
    input:
        fastq1 = "fastq/DNA/{sample}_R1.fastq.gz"
        fastq2 = "fastq/DNA/{sample}_R2.fastq.gz"
    output:
        bam = "DNA_bam/{sample}-ready.bam"
    params:
        config["reference"]["ref"],
        known_indels = "/data/ref_genomes/hg19/variation/Mills_and_1000G_gold_standard.indels.vcf.gz",
        out_recal = "bam/{sample}.recal_gpu.txt"
    log:
        "logs/GPU_align/{sample}.log"
    threads:
        2
    singularity:
        config["singularitys"]["GPU"]
    shell:
        "(pbrun fq2bam "
        "--ref {params.ref} "
        "--in-fq {input.fastq1} {input.fastq2} "
        "--knownSites {params.known_indels} "
        "--out-bam {output.bam} "
        "--out-recal-file {params.out_real} "
        "-tmp-dir GPU_run "
        ") &> {log}"


pbrun fq2bam --ref Ref/Homo_sapiens_assembly38.fasta \
--in-fq S1_1.fastq.gz S1_2.fastq.gz \
--knownSites Ref/Homo_sapiens_assembly38.known_indels.vcf.gz \
--out-bam mark_dups_gpu.bam \
--out-recal-file recal_gpu.txt \
--tmp-dir /raid/myrun
