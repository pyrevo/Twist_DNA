import math
from subprocess import check_output


in_fastq_ref = "/data/ref_genomes/hg19/bwa/BWA_0.7.10_refseq/hg19.with.mt.fasta"
# in_vcf = open("/home/jonas/Investigations/Twist_DNA/BMS53sonic-ensemble.final.no.introns.AD20.ffpe.tsv")
# out_vcf = open("/home/jonas/Investigations/Twist_DNA/BMS53sonic-ensemble.final.no.introns.AD20.ffpe.Mbp.tsv", "w")
in_vcf = open("/scratch/wp1/nobackup/ngs/utveckling/analys/2021/20210112_HN_Kapa_Twist_GMS560v1/Results/DNA/BMS53sonic/vcf/\
BMS53sonic-ensemble.final.no.introns.AD20.ffpe.tsv")
out_vcf = open("/scratch/wp1/nobackup/ngs/utveckling/analys/2021/20210112_HN_Kapa_Twist_GMS560v1/Results/DNA/BMS53sonic/vcf/\
BMS53sonic-ensemble.final.no.introns.AD20.ffpe.Mbp.tsv", "w")


AA_dict = {"Ala": ["GCT", "GCC", "GCA", "GCG"],
           "Ile": ["ATT", "ATC", "ATA"],
           "Arg": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
           "Leu": ["CTT", "CTC", "CTA", "CTG", "TTA", "TTG"],
           "Asn": ["AAT", "AAC"],
           "Lys": ["AAA", "AAG"],
           "Asp": ["GAT", "GAC"],
           "Met": ["ATG"],
           "Phe": ["TTT", "TTC"],
           "Cys": ["TGT", "TGC"],
           "Pro": ["CCT", "CCC", "CCA", "CCG"],
           "Gln": ["CAA", "CAG"],
           "Ser": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
           "Glu": ["GAA", "GAG"],
           "Thr": ["ACT", "ACC", "ACA", "ACG"],
           "Trp": ["TGG"],
           "Gly": ["GGT", "GGC", "GGA", "GGG"],
           "Tyr": ["TAT", "TAC"],
           "His": ["CAT", "CAC"],
           "Val": ["GTT", "GTC", "GTA", "GTG"],
           "*":	["TAA", "TGA", "TAG"]
           }

dna_opposite = {"A": "T", "T": "A", "G": "C", "C": "G", "X": "X"}


header = True
candidate_list = []
prev_pos = 0
prev_chrom = ""
prev_lline = []
for line in in_vcf:
    if header:
        if line[:6] == "#CHROM":
            header = False
        out_vcf.write(line)
        continue
    lline = line.strip().split()
    chrom = lline[0]
    pos = int(lline[1])
    ref = lline[3]
    alt = lline[4]
    if len(ref) > 1 or len(alt) > 1:
        continue
    if chrom == prev_chrom and pos - prev_pos <= 2:
        if not (candidate_list != [] and candidate_list[-1][0] == prev_chrom and int(candidate_list[-1][1]) == prev_pos):
            candidate_list.append(prev_lline)
        candidate_list.append(lline)
    prev_pos = pos
    prev_chrom = chrom
    prev_lline = lline

prev_aa_nr = 0
prev_chrom = ""
prev_pos = 0
prev_candidate = []
Multibp_list = []
for candidate in candidate_list:
    gene_change = candidate[7].split("|c.")[1].split("|")[0]
    gene_pos = int(gene_change[:-3])
    aa_nr = math.ceil(gene_pos / 3.0)
    print(candidate[:2], gene_change, aa_nr)
    chrom = candidate[0]
    pos = int(candidate[1])
    if chrom == prev_chrom and pos - prev_pos <= 2 and aa_nr == prev_aa_nr:
        if Multibp_list != []:
            print(Multibp_list[-1][1][0:2], chrom, pos, prev_chrom, prev_pos)
        if Multibp_list != [] and Multibp_list[-1][1][0] == prev_chrom and int(Multibp_list[-1][1][1]) == prev_pos:
            Multibp_list[-1].append(candidate)
        else:
            Multibp_list.append([prev_candidate, candidate])
    prev_aa_nr = aa_nr
    prev_chrom = chrom
    prev_pos = pos
    prev_candidate = candidate

for Multibp in Multibp_list:
    codon_pos_list = []
    for variant in Multibp:
        gene_change = variant[7].split("|c.")[1].split("|")[0]
        gene_pos = int(gene_change[:-3])
        codon_pos = (gene_pos - 1) % 3
        codon_pos_list.append(codon_pos)
    flip_bp = False
    if codon_pos_list[0] > codon_pos_list[1]:
        flip_bp = True

    ref = ["X", "X", "X"]
    alt = ["X", "X", "X"]
    for variant in Multibp:
        gene_change = variant[7].split("|c.")[1].split("|")[0]
        dna_change = gene_change[-3:].split(">")
        gene_pos = int(gene_change[:-3])
        codon_pos = (gene_pos - 1) % 3
        if flip_bp:
            ref[codon_pos] = dna_opposite[dna_change[0]]
            alt[codon_pos] = dna_opposite[dna_change[1]]
        else:
            ref[codon_pos] = dna_change[0]
            alt[codon_pos] = dna_change[1]

    print(codon_pos_list)
    chrom = Multibp[0][0]
    if flip_bp:
        pos = int(Multibp[0][1]) - (2 - codon_pos_list[0])
    else:
        pos = int(Multibp[0][1]) - codon_pos_list[0]
    i = 0
    for bp in ref:
        if bp == "X":
            fasta_pos = pos + i
            command = "samtools faidx " + in_fastq_ref + " " + chrom + ":" + str(fasta_pos) + "-" + str(fasta_pos)
            ref_bp = check_output(command, shell=True).decode("utf-8").split("\n")[1]
            ref[i] = ref_bp
            alt[i] = ref_bp
        i += 1
    out_vcf.write(chrom + "\t" + str(pos) + "\t.\t" + "".join(ref) + "\t" + "".join(alt) + "\t.\tPASS\n")


# print()
# for Multibp in Multibp_list:
#     codon_dict = {0:"X", 1:"X", 2:"X"}
#     wt_codon_dict = {0:"X", 1:"X", 2:"X"}
#     codon_pos_list = []
#     for variant in Multibp:
#         print(variant)
#         gene_change = variant[7].split("|c.")[1].split("|")[0]
#         dna_change = gene_change[-3:].split(">")
#         gene_pos = int(gene_change[:-3])
#         aa_nr = math.ceil(gene_pos / 3.0)
#         codon_pos = (gene_pos - 1) % 3
#         protein_change = variant[7].split("|p.")[1].split("|")[0]
#         wt_aa = protein_change[:3]
#         print(variant[:2], gene_change, dna_change, gene_pos, aa_nr, protein_change, wt_aa, codon_pos)
#         codon_pos_list.append(codon_pos)
#         codon_dict[codon_pos] = dna_change[1]
#         wt_codon_dict[codon_pos] = dna_change[0]
#     if codon_pos_list[0] > codon_pos_list[1]:
#         for codon_pos in codon_dict:
#             codon_dict[codon_pos] = dna_opposite[codon_dict[codon_pos]]
#     print(codon_dict, wt_codon_dict)
#     print(AA_dict[wt_aa])
#     wt_codon = ""
#     i = 0
#     while i <= 2:
#         wt_codon += wt_codon_dict[i]
#         i += 1
#     best_codon = ""
#     for codon in AA_dict[wt_aa]:
#         print(codon, wt_codon)
#         i = 0
#         match = True
#         for base in codon:
#             if not (base == wt_codon[i] or wt_codon[i] == "X"):
#                 match = False
#             i += 1
#         if match:
#             if best_codon != "":
#                 print("Ops!")
#                 quit()
#             best_codon = codon
#     wt_codon = best_codon
#     print(wt_codon)
