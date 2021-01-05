
rule freebayes:
    input:
        bam = "DNA_bam/{sample}-ready.bam",
        bai = "DNA_bam/{sample}-ready.bam.bai",
        ref = config["reference"]["ref"],
        bed = config["bed"]["bedfile"]
    output:
        temp("freebayes/{sample}.freebayes.unsort.vcf")
    log:
        "logs/variantCalling/freebayes/{sample}.log"
    params:
        freebayes_singularity = config["singularity"]["execute"] + config["singularity"]["freebayes"],
        bcftools_singularity = config["singularity"]["execute"] + config["singularity"]["bcftools"],
        extra = " --min-alternate-fraction 0.01 --allele-balance-priors-off --pooled-discrete --pooled-continuous --report-genotype-likelihood-max --genotype-qualities --strict-vcf --no-partial-observations "
    shell:
        "({params.freebayes_singularity} freebayes {params.extra} -t {input.bed} -f {input.ref} {input.bam} |"
        " {params.bcftools_singularity} bcftools filter -i 'ALT=\"<*>\" || QUAL > 5' |"
        " awk -F$'\t' -v OFS='\t' '{{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDXkmryswbvhdx]/, \"N\", $4) }} {{print}}' > {output}) &> {log}"


rule sortFreebayes:
    input:
        "freebayes/{sample}.freebayes.unsort.vcf"
    output:
        temp("freebayes/{sample}.freebayes.fixAF.vcf")
    singularity:
        config["singularity"]["bcftools"]
    log:
        "logs/variantCalling/freebayes/{sample}.sort.log"
    shell:
        "(bcftools sort -o {output} -O v {input}) &> {log}"
