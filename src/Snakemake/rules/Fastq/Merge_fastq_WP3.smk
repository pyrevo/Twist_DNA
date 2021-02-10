
L_numbers = ["1", "2", "3", "4"]
S_dna = []
for s in config["DNA_Samples"].values():
    S_dna.append(s)
fastq1_files = ["fastq_temp/" + s + "_" + i + "_L00" + L + "_R1_001.fastq.gz" for s, i, L in zip(config["DNA_Samples"], S_dna, L)]
fastq2_files = ["fastq_temp/" + s + "_" + i + "_L00" + L + "_R2_001.fastq.gz" for s, i, L in zip(config["DNA_Samples"], S_dna, L)]

rule merge_Fastq:
    input:
        fastq1_files = fastq1_files,
        fastq2_files = fastq2_files,
    output:
        fastq1 = "fastq/{sample}_R1.fastq.gz",
        fastq2 = "fastq/{sample}_R2.fastq.gz",
    log:
        "logs/fastq/merge/{sample}.log",
    shell:
        "zcat {input.fastq1_files} | pigz > {output.fastq1} && "
        "zcat {input.fastq2_files} | pigz > {output.fastq2}"
