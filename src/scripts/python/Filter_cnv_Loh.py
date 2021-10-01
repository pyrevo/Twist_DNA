
import glob
import gzip
import os
import subprocess
import logging

if len(snakemake.log) > 0:
    logging.basicConfig(filename=snakemake.log[0], level=logging.INFO)

cnv_purity = snakemake.params.purity
cnv_relevant_genes = open(snakemake.input.relevant_genes)
cnvkit_files = snakemake.input.cnvkit_segments
GATK_CNV_files = snakemake.input.GATK_CNV_segments
cnv_bed_file = open(snakemake.input.bed_file)
cnv_relevant = open(snakemake.output.relevant_cnvs, "w")
cnv_relevant_clinical = open(snakemake.output.relevant_cnvs_clinical, "w")
in_path = snakemake.params.in_path
out_path = snakemake.params.out_path

cnv_relevant_clinical.write("Method\tsample\tgenes\tchrom\tstart_pos\tend_pos\tCN_CNVkit\tCN_GATK\tGATK_BAF\tregion_size")
cnv_relevant_clinical.write("\tnr_exons\tnr_obs\tdepth\tpurity\n")
cnv_relevant.write("Method\tsample\tgenes\tchrom\tstart_pos\tend_pos\tCN\tBAF\tregion_size\tnr_exons\tnr_obs_cov")
cnv_relevant.write("\tnr_obs_baf\tpurity\n")

chrom_len = {"chr1": 249250621, "chr2": 243199373, "chr3": 198022430, "chr4": 191154276, "chr5": 180915260, "chr6": 171115067,
             "chr7": 159138663, "chr8": 146364022, "chr9": 141213431, "chr10": 135534747, "chr11": 135006516, "chr12": 133851895,
             "chr13": 115169878, "chr14": 107349540, "chr15": 102531392, "chr16": 90354753, "chr17": 81195210, "chr18": 78077248,
             "chr19": 59128983, "chr20": 63025520, "chr21": 48129895, "chr22": 51304566, "chrX": 155270560, "chrY": 59373566}

relevant_genes = {}
for line in cnv_relevant_genes:
    relevant_genes[line.strip()] = {}

gene_regions = {}
for line in cnv_bed_file:
    lline = line.strip().split("\t")
    chrom = lline[0]
    start = lline[1]
    end = lline[2]
    name = lline[3]
    gene = name.split("_")[0]
    if gene in gene_regions:
        gene_regions[gene][2] = end
    else:
        gene_regions[gene] = [chrom, start, end]


'''Pathological purity'''
sample_purity_dict = {}
for row in cnv_purity:
    column = row.split(";")
    purity = float(column[1])
    if purity == 0:
        purity = 1.0
    sample_purity_dict[column[0]] = [0, 0, 0, purity]


'''Extract events from CNVkit'''
CNVkit_regions = []
for cnv_file_name in cnvkit_files:
    seg = open(cnv_file_name)
    sample = cnv_file_name.split("/")[-1].split(".")[0]
    sample2 = sample.split("-LoH")[0]
    next(seg)
    for line in seg:
        lline = line.strip().split("\t")
        if lline[3].find("MUC6") != -1:
            continue
        CNVkit_CR = float(lline[4])
        CNVkit_MAF = lline[5]
        if CNVkit_MAF == "":
            CNVkit_MAF = 0.0
        else:
            CNVkit_MAF = float(CNVkit_MAF)
        depth = "0"
        nr_probes = "0"
        if len(lline) >= 12 :
            depth = lline[11]
            nr_probes = lline[12]
        length = int(lline[2]) - int(lline[1]) + 1
        purity = sample_purity_dict[sample2][3]
        CN_CNVkit = round(2*pow(2, CNVkit_CR), 2)
        if CN_CNVkit < 1.2 and abs(CNVkit_MAF - 0.5) > 0.15 and length > 100000:
            CNVkit_regions.append([
                "CNVkit", sample2, lline[0], lline[1], lline[2], CNVkit_CR, CNVkit_MAF, purity, length, CN_CNVkit, lline[3],
                depth, nr_probes
            ])
        elif CN_CNVkit > 4.0 and length > 10000:
            CNVkit_regions.append([
                "CNVkit", sample2, lline[0], lline[1], lline[2], CNVkit_CR, CNVkit_MAF, purity, length, CN_CNVkit, lline[3],
                depth, nr_probes
            ])
    seg.close()

'''Extract events from GATK'''
GATK_regions = []
for cnv_file_name in GATK_CNV_files:
    seg = open(cnv_file_name)
    sample = cnv_file_name.split("/")[-1].split(".")[0]
    sample2 = sample.split("_clean")[0]
    header = True
    for line in seg:
        if header:
            if line[:6] == "CONTIG":
                header = False
            continue
        lline = line.strip().split("\t")
        GATK_CR = float(lline[6])
        GATK_MAF = lline[8]
        if GATK_MAF == "NaN":
            GATK_MAF = 0.5
        else:
            GATK_MAF = float(GATK_MAF)
        length = int(lline[2]) - int(lline[1]) + 1
        Points_CR = int(lline[3])
        Points_MAF = int(lline[4])
        CN_GATK_100 = round(2*pow(2, GATK_CR), 2)
        purity = sample_purity_dict[sample2][3]
        GATK_corrected_CN = round(2 + (CN_GATK_100 - 2) * (1/purity), 1)
        if GATK_corrected_CN < 1.2 and abs(GATK_MAF - 0.5) > 0.15 and length > 100000 and (Points_CR > 20 or Points_MAF > 20):
            GATK_regions.append([
                "GATK_CNV", sample2, lline[0], lline[1], lline[2], GATK_CR, GATK_MAF, purity, length,
                GATK_corrected_CN, str(Points_CR), str(Points_MAF)
            ])
        elif GATK_corrected_CN > 4.0 and length > 10000:
            GATK_regions.append([
                "GATK_CNV", sample2, lline[0], lline[1], lline[2], GATK_CR, GATK_MAF, purity, length,
                GATK_corrected_CN, str(Points_CR), str(Points_MAF)
            ])
    seg.close()

'''Add GATK CN and BAF to CNVkit clinical results'''
Both_regions = []
for region1 in GATK_regions:
    found = False
    for region2 in CNVkit_regions:
        if region1[1] == region2[1] and region1[2] == region2[2]:
            if (
                (int(region1[3]) >= int(region2[3]) and int(region1[3]) <= int(region2[4])) or
                (int(region1[4]) >= int(region2[3]) and int(region1[4]) <= int(region2[4])) or
                (int(region1[3]) <= int(region2[3]) and int(region1[4]) >= int(region2[4]))
            ):
                if (region1[9] < 1.5 and region2[9] < 1.5) or (region1[9] > 2.5 and region2[9] > 2.5):
                    if len(Both_regions) > 0:
                        if not (region2[1] == Both_regions[-1][1] and region2[2] == Both_regions[-1][2] and
                                region2[3] == Both_regions[-1][3] and region2[4] == Both_regions[-1][4]):
                            Both_regions.append(region2 + [region1[9]] + [region1[6]])
                    else:
                        Both_regions.append(region2 + [region1[9]] + [region1[6]])
                    found = True
    if found:
        Both_regions.append(region1)

for region in Both_regions:
    method = region[0]
    genes = ""
    nr_exons = ""
    gene_dict = {}
    if method == "CNVkit":
        exons = region[10]
        nr_exons = len(exons.split(","))
        for exon in exons.split(","):
            gene = exon.split("_")[0]
            gene_dict[gene] = ""
        for gene in gene_dict:
            if genes == "":
                genes = gene
            else:
                genes += "," + gene
    nr_obs_cov = ""
    nr_obs_baf = ""
    if method == "GATK_CNV":
        nr_obs_cov = region[10]
        nr_obs_baf = region[11]
    cnv_relevant.write(method + "\t" + region[1] + "\t" + genes + "\t" + region[2] + "\t" + region[3] + "\t" +
                       region[4] + "\t" + str(region[9]) + "\t" + str(region[6]) + "\t" + str(region[8]) + "\t" +
                       str(nr_exons) + "\t" + nr_obs_cov + "\t" + nr_obs_baf + "\t" + str(region[7]) + "\n")
    if method == "CNVkit":
        found_gene = False
        clinical_gene = ""
        for gene in genes.split(","):
            if gene in relevant_genes:
                found_gene = True
                clinical_gene = gene
        '''change this'''
        if (region[9] < 1.5 and (region[2] == "chr1" or region[2] == "chr19")):
            clinical_gene = "1p19q?"
            cnv_relevant_clinical.write(
                method + "\t" + region[1] + "\t" + clinical_gene + "\t" + region[2] + "\t" + region[3] + "\t" +
                region[4] + "\t" + str(region[9]) + "\t" + str(region[13]) + "\t" + str(region[14]) +
                "\t" + str(region[8]) + "\t" + str(nr_exons) + "\t" + region[12] + "\t" + region[11] + "\t" + str(region[7]) + "\n"
            )
        if (found_gene and region[9] > 2.5) or (region[9] < 1.5 and region[2] == "chr1"):
            cnv_relevant_clinical.write(
                method + "\t" + region[1] + "\t" + clinical_gene + "\t" + region[2] + "\t" + region[3] + "\t" +
                region[4] + "\t" + str(region[9]) + "\t" + str(region[13]) + "\t" + str(region[14]) +
                "\t" + str(region[8]) + "\t" + str(nr_exons) + "\t" + region[12] + "\t" + region[11] + "\t" + str(region[7]) + "\n"
            )
