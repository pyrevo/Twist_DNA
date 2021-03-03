# rule fastqc_bam:
#     input:
#         #"fastq/{sample}_R1-fastq.gz" ##one for each R1 and one for R2 should be from a samples.yaml file
#         #"DNA_bam/{sample}-ready.bam"
#         #"bam/{sample}-sort.bam"
#         bam = "STAR2/{sample}Aligned.sortedByCoord.out.bam"
#     output:
#         html="qc/{sample}/{sample}Aligned.sortedByCoord.out_fastqc.html",
#         zip="qc/{sample}/{sample}Aligned.sortedByCoord.out_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
#     params:
#         outdir = "qc/{sample}/"
#     log:
#         "logs/qc/fastqc/{sample}Aligned.sortedByCoord.log"
#     threads: 10
#     singularity:
#         config["singularity"]["fastqc"]
#     shell:
#         "(fastqc --quiet -t {threads} --outdir {params.outdir} {input}) &> {log}"
#     # wrapper:
#     #     "0.38.0/bio/fastqc"


rule fastqcR1:
    input:
        "fastq/DNA/{sample}_R1.fastq.gz",  ##one for each R1 and one for R2 should be from a samples.yaml file
    output:
        html="qc/{sample}/{sample}_R1_fastqc.html",
        zip="qc/{sample}/{sample}_R1_fastqc.zip",  # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        outdir="qc/{sample}/",
    log:
        "logs/qc/fastqc/{sample}_R1.log",
    singularity:
        config["singularity"].get("fastqc", config["singularity"].get("default", ""))
    shell:
        "(fastqc --quiet --outdir {params.outdir} {input}) &> {log}"


rule fastqcR2:
    input:
        "fastq/DNA/{sample}_R2.fastq.gz",  ##one for each R1 and one for R2 should be from a samples.yaml file
    output:
        html="qc/{sample}/{sample}_R2_fastqc.html",
        zip="qc/{sample}/{sample}_R2_fastqc.zip",  # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        outdir="qc/{sample}/",
    log:
        "logs/qc/fastqc/{sample}_R2.log",
    singularity:
        config["singularity"].get("fastqc", config["singularity"].get("default", ""))
    shell:
        "(fastqc --quiet --outdir {params.outdir} {input}) &> {log}"
