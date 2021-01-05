
rule vep:
    input:
        vcf = "recall/{sample}.ensemble.vcf.gz",
        cache = config["configfiles"]["vep"], #"/opt/vep/.vep", ## always remeber the --bind vep-data:/opt/vep/.vep command in singularity args
        fasta = config["reference"]["ref"]
    output:
        vcf = temp("variantCalls/annotation/raw/{sample}_{seqID}.raw.vcf")
    params:
        "--check_existing --pick --sift b --polyphen b --ccds --uniprot --hgvs --symbol --numbers --domains --regulatory --canonical --protein --biotype --uniprot --tsl --appris --gene_phenotype --af --af_1kg --af_gnomad --max_af --pubmed --variant_class "
        # "--everything --check_existing --pick"  #--exclude_null_alleles
    log:
        "logs/variantCalling/vep/{sample}.log"
    singularity:
        config["singularitys"]["vep"]
    threads: 10
    shell:
        "(vep --vcf --no_stats -o {output.vcf} -i {input.vcf} --dir_cache {input.cache} --fork {threads} --cache --refseq --offline --fasta {input.fasta} {params} ) &> {log}"
