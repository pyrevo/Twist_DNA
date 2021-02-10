
L_numbers = ["1", "2", "3", "4"]
S_numbers = ["1", "2", "3", "4", "5", "6", "7", "8"]

rule merge_Fastq:
    input:
        fastq1_files = expand("fastq_temp/{{sample}}_S{S_number}_L00{L_number}_R1_001.fastq.gz", L_number=L_numbers, S_number=S_numbers),
        fastq2_files = expand("fastq_temp/{{sample}}_S{S_number}_L00{L_number}_R2_001.fastq.gz", L_number=L_numbers, S_number=S_numbers),
    output:
        fastq1 = "fastq/{sample}_R1.fastq.gz",
        fastq2 = "fastq/{sample}_R2.fastq.gz",
    log:
        "logs/fastq/merge/{sample}.log",
    shell:
        "zcat {input.fastq1_files} | pigz > {output.fastq1} && "
        "zcat {input.fastq2_files} | pigz > {output.fastq2}"
