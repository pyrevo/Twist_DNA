
rule varscan:
    input:
        bam="Bam/DNA/{sample}-ready.bam",
        bai="Bam/DNA/{sample}-ready.bam.bai",
        ref=config["reference"]["ref"],
        bed=config["bed"]["bedfile"],
    output:
        temp("varscan/{sample}.varscan.vcf"),
    params:
        samtools_singularity=config["singularity"]["execute"] + config["singularity"].get(
            "samtools", config["singularity"].get("default", "")
        ),
        varscan_singularity=config["singularity"]["execute"] + config["singularity"].get("varscan", config["singularity"].get("default", "")),
        mpileup="-d 1000 -L 1000",
        varscan="--min-coverage 5 --p-value 0.98 --strand-filter 1 --min-var-freq 0.01 --output-vcf --variants",
    log:
        "logs/variantCalling/varscan/{sample}.log",
    shell:
        # Exchange alleles frequences from ex 10.5% to 0.105
        # Exchange alleles frequences from ex 10% to 0.10
        # Exchange alleles frequences from 100% to 1.00
        # " sed 's/FREQ/AF/' " #Exchange FORMAT/FREQ to FORMAT/AF
        "({params.samtools_singularity} samtools mpileup -f {input.ref} {params.mpileup} -l {input.bed} {input.bam} |"
        " grep -v -P '\t0\t\t$' |"
        " {params.varscan_singularity} java -jar /usr/local/share/varscan-2.4.3-0/VarScan.jar mpileup2cns {params.varscan} |"
        " awk -F$'\t' -v OFS='\t' '{{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDXkmryswbvhdx]/, \"N\", $4) }} {{print}}' | "
        " awk -F$'\t' -v OFS='\t' '{{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDXkmryswbvhdx]/, \"N\", $5) }} {{print}}' | "
        " sed 's/:\([0-9]\+\)\.\([0-9]\+\)%:/:0.\\1\\2:/' | "
        " sed 's/:\([0-9]\+\)%:/:0.\\1:/' | "
        " sed 's/:100%:/:1.00:/' "
        " > {output}) 2> {log}"


rule sortVarscan:
    input:
        "varscan/{sample}.varscan.vcf",
    output:
        temp("varscan/{sample}.varscan.fixAF.vcf"),
    singularity:
        config["singularity"].get("bcftools", config["singularity"].get("default", ""))
    log:
        "logs/variantCalling/varscan/{sample}.sort.log",
    shell:
        "(bgzip {input} && tabix {input}.gz && "
        "bcftools sort -o {output} -O v {input}.gz) &> {log}"
