
import re

fusions = open(snakemake.input.fusions)
genes = open(snakemake.input.genes)
fusions_filtered = open(snakemake.output.fusions, "w")


druggable_gene_dict = {}
for line in genes :
    if line[0] == ">" :
        gene = line[1:].split("_")[0]
        druggable_gene_dict[gene] = ""

header = True
for line in fusions :
    if header :
        fusions_filtered.write(line)
        header = False
        continue
    column = line.strip().split("\t")
    gene1 = column[11]
    gene2 = column[13]
    gene1_list = re.split('_|&', gene1)
    gene2_list = re.split('_|&', gene2)
    gene_list = gene1_list + gene2_list
    found_gene = False
    for gene in gene_list :
        if gene in druggable_gene_dict :
            found_gene = True
    if found_gene :
        fusions_filtered.write(line)
