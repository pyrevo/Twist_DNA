
import subprocess
import sys

bed = snakemake.input.bed
bam_file = snakemake.input.bam
outfile = open(snakemake.output.coverage, "w")
outfile2 = open(snakemake.output.coverage2, "w")

header = "#Chr\tStart_hg19\tEnd_hg19\tGene\tCDS_mut_syntax\tAA_mut_syntax\tReport"
header += "\tcomment\tExon\tAccession_number\tCoverage\tPosition\n"

outfile.write(header)
outfile2.write(header)


'''Find positions to report and gene regions to analyse'''
inv_pos = {}
gene_regions = []
header = True
prev_gene = ""
prev_chrom = ""
first_gene = True
gene_start_pos = 0
gene_end_pos = 0
for line in bed:
    if header:
        header = False
        continue
    lline = line.strip().split("\t")
    Report = lline[6]
    if Report == "region":
        continue
    chrom = lline[0].split(".")[0].split("0")[-1]
    start_pos = lline[1]
    end_pos = lline[2]
    gene = lline[3]
    pos = int(start_pos)
    while pos <= int(end_pos):
        inv_pos[chrom + "_" + str(pos)] = lline
        pos += 1
    if first_gene:
        prev_gene = gene
        prev_chrom = chrom
        gene_start_pos = start_pos
        gene_end_pos = end_pos
        first_gene = False
        continue
    if prev_gene != gene:
        gene_regions.append([int(prev_chrom), "chr" + prev_chrom + ":" + gene_start_pos + "-" + gene_end_pos])
        prev_gene = gene
        prev_chrom = chrom
        gene_start_pos = start_pos
        gene_end_pos = end_pos
    else:
        gene_end_pos = end_pos
gene_regions.append([int(chrom), "chr" + chrom + ":" + gene_start_pos + "-" + gene_end_pos])

'''Remove identical regions'''
gene_regions.sort()
gene_regions_temp = []
for gene_region in gene_regions:
    gene_region = gene_region[1]
    if gene_region not in gene_regions_temp:
        gene_regions_temp.append(gene_region)
gene_regions = gene_regions_temp

'''Report all interesting positions (All but region) with coverage < 50'''
depth_dict = {}
for region in gene_regions:
    sample = bam_file.split("/")[-1].split(".")[0]
    outfile_name = "DATA/gene_depth_" + sample + ".txt"
    subprocess.call("samtools depth -a -r " + region + " " + bam_file + " > " + outfile_name, shell=True)
    depth_file = open(outfile_name)
    for line in depth_file:
        lline = line.strip().split("\t")
        chrom = lline[0][3:]
        pos = lline[1]
        key = chrom + "_" + pos
        if key in inv_pos:
            coverage = int(lline[2])
            if coverage < 50:
                for info in inv_pos[key]:
                    outfile.write(info + "\t")
                outfile.write(str(coverage) + "\t" + pos + "\n")
            for info in inv_pos[key]:
                outfile2.write(info + "\t")
            outfile2.write(str(coverage) + "\t" + pos + "\n")
    depth_file.close()
outfile.close()
outfile2.close()
