
#snakemake -p -j 64 --drmaa "-A wp1 -p core -n {cluster.n} -t {cluster.time}"  -s ./Twist_exome_WP3.smk --use-singularity --singularity-args "--bind /data --bind /projects --bind /scratch " --cluster-config Config/Slurm/cluster.json

configfile: "Twist_exome_WP3.yaml"

wildcard_constraints:
    unit="[A-Za-z0-9-]+",
    sample="[^.]+",

def get_input():
    input_list = []
    '''Fastq'''
    input_list.append(["fastq/" + s + "_R1.fastq.gz" for s in config["DNA_Samples"]])
    input_list.append(["fastq/" + s + "_R2.fastq.gz" for s in config["DNA_Samples"]])

    '''Alignment'''
    input_list.append(["alignment/" + s + ".bam" for s in config["DNA_Samples"]])
    input_list.append(["alignment/" + s + ".bam.bai" for s in config["DNA_Samples"]])

    '''SNV / INDEL'''
    input_list.append(["haplotypecaller/" + s + ".vcf.gz" for s in config["DNA_Samples"]])
    input_list.append(["haplotypecaller/" + s + ".vep.vcf.gz" for s in config["DNA_Samples"]])
    input_list.append(["haplotypecaller/" + s + ".vep.filteredSNP.filteredINDEL.vcf.gz" for s in config["DNA_Samples"]])

    return input_list

rule all:
    input:
        get_input()


include: "src/Snakemake/workflow/Twist_exome_WP3_workflow.smk"
