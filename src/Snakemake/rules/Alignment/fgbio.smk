

import src.lib.python.utils as utils

def get_now():
    from datetime import datetime

    return datetime.now().strftime("%Y%m%d")

_bwa_mem_fgbio1_input = ["fastq/DNA/{sample}_R1.fastq.gz", "fastq/DNA/{sample}_R2.fastq.gz"]
_temp_bwa_mem_fgbio1_output = "alignment/{sample}.prep_fgbio.sort.bam"
_bwa_men_fgbio1_log = "logs/map/bwa/{sample}.log"
_bwa_benchmark_fgbio1 = "benchmarks/bwa/mem/{sample}.tsv"
_pu = "{sample}"
_prep_fgbio1_input = "alignment/{sample}.fgbio.sort.bam"
if "units" in config:
    _bwa_mem_fgbio1_input = ["fastq/DNA/{sample}_{unit}_R1.fastq.gz", "fastq/DNA/{sample}_{unit}_R2.fastq.gz"]
    _temp_bwa_mem_fgbio1_output = "alignment/{sample}_{unit}.prep_fgbio.sort.bam"
    _bwa_men_fgbio1_log = "logs/map/bwa/{sample}_{unit}.log"
    _bwa_mem_fgbio1_benchmark = "benchmarks/bwa/mem/{sample}_{units}.tsv"
    _pu = "{sample}_{unit}"
    _prep_fgbio1_input = "alignment/{sample}.merged.prep_fgbio.sort.bam"

try:
    _bwa_mem_fgbio1_input = bwa_mem_fgbio1_input
except:
    pass

rule bwa_mem_fgbio:
    input:
        reads=_bwa_mem_fgbio1_input,
    output:
        bam=temp(_temp_bwa_mem_fgbio1_output),
    log:
        _bwa_men_fgbio1_log,
    params:
        index=config["reference"]["ref"],
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:illumina\tPU:" + _pu + ' -v 1 ' + config.get("bam_extra", ""),
        sort="samtools",
        sort_order="coordinate",
        sort_extra="-@ 10",
    threads: 10
    #benchmark:
    #    repeat(_bwa_benchmark, config.get("benchmark", {}).get("repeats", 1))
    singularity:
        config["singularity"].get("bwa", config["singularity"].get("default", ""))
    wrapper:
        "0.70.0/bio/bwa/mem"


rule merge_bam_for_fgbio:
    input:
        lambda wildcards: ["alignment/" + wildcards.sample + "_" + unit+ ".prep_fgbio.sort.bam" for unit in utils.get_units(units, wildcards.sample)]
    output:
        temp("alignment/{sample}.merged.prep_fgbio.sort.bam"),
    singularity:
        config["singularity"].get("samtools", config["singularity"].get("default", ""))
    shell:
        """
        samtools merge -c -p {output} {input}
        """


rule prep_fgbio1:
    input:
        bam=_prep_fgbio1_input,
    output:
        temp("alignment/{sample}.prep_fgbio.sort.bam"),
    log:
        "logs/fgbio/bwa1/{sample}.log",
    params:
        bwa_singularity=config["singularity"]["execute"] + config["singularity"].get(
            "bwa", config["singularity"].get("default", "")
        ),
        bamsormadup_singularity=config["singularity"]["execute"] + config["singularity"]["bamsormadup"],
        umis_singularity=config["singularity"]["execute"] + config["singularity"]["umis"],
        samtools_singularity=config["singularity"]["execute"] + config["singularity"].get(
            "samtools", config["singularity"].get("default", "")
        ),
        index=config["reference"]["ref"],
        extra=r"-c 250 -M -R '@RG\tID:{sample}\tSM:{sample}\tPL:illumina\tPU:{sample}' -v 1",
        tmp_dir="tmpfile=bam/{sample}",
    threads: 10
    shell:
        "({params.bamsormadup_singularity} bamsormadup {params.tmp_dir} inputformat=bam threads={threads} outputformat=bam level=0 SO=coordinate {input.bam}"
        " | {params.umis_singularity} umis bamtag -"
        " | {params.samtools_singularity} samtools view -b -o {output} - ) &> {log}"


# rule samtools_index_fgbio1:
#     input:
#         "bam/{sample}-sort.bam",
#     output:
#         "bam/{sample}-sort.bam.bai",
#     log:
#         "logs/fgbio/samtools_index1/{sample}.log",
#     singularity:
#         config["singularity"]["samtools"]
#     shell:
#         "(samtools index {input} {output}) &> {log}"


rule fgbio:
    input:
        bam="alignment/{sample}.prep_fgbio.sort.bam",
        ref=config["reference"]["ref"],
    output:
        fq1=temp("fastq_temp/{sample}-cumi-R1.fq.gz"),
        fq2=temp("fastq_temp/{sample}-cumi-R2.fq.gz"),
        qc="qc/{sample}/{sample}_fgbio.txt",
    params:
        bam_tmp="alignment/{sample}-cumi-1-bamtofastq-tmp",
        fgbio_singularity=config["singularity"]["execute"] + config["singularity"]["fgbio"],
        bamtofastq_singularity=config["singularity"]["execute"] + config["singularity"]["bamtofastq"],
    log:
        "logs/fgbio/{sample}.log",
    shell:
        "{params.fgbio_singularity} fgbio GroupReadsByUmi -i {input.bam} --edits=1 --min-map-q=1 -t RX -s adjacency -f {output.qc}"
        " | {params.fgbio_singularity} fgbio CallMolecularConsensusReads -i /dev/stdin -o /dev/stdout"
        " --min-input-base-quality=2 --min-reads=1 --max-reads=100000 --output-per-base-tags=false --sort-order=:none:"
        " | {params.fgbio_singularity} fgbio FilterConsensusReads -i /dev/stdin -o /dev/stdout -r {input.ref} --min-reads=1 --min-base-quality=13"
        " | {params.bamtofastq_singularity} bamtofastq collate=1 tags=cD,cM,cE gz=1 T={params.bam_tmp} F={output.fq1} F2={output.fq2}"
        " > {log}"


rule bwa_mem_fgbio2:
    input:
        reads=["fastq_temp/{sample}-cumi-R1.fq.gz", "fastq_temp/{sample}-cumi-R2.fq.gz"],
    output:
        "Bam/DNA/{sample}-ready.bam",
    log:
        "logs/fgbio/bwa2/{sample}.log",
    params:
        index=config["reference"]["ref"],
        extra=r"-C -c 250 -M -R '@RG\tID:{sample}\tSM:{sample}\tPL:illumina\tPU:{sample}' -v 1",
    threads: 10
    singularity:
        config["singularity"].get("bwa", config["singularity"].get("default", ""))
    shell:
        "(bwa mem -t {threads} {params.extra} {params.index} {input.reads} | samtools sort -@ {threads} -m 3G -o {output} - ) &> {log}"
