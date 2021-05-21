
vcf = open(snakemake.input.vcf)
artifacts = open(snakemake.input.artifacts)
output_tmb = open(snakemake.input.tmb, "w")


FFPE_SNV_artifacts = {}
header = True
for line in artifacts:
    if header:
        continue
    lline = line.strip().split("\t")
    chrom = lline[0]
    pos = lline[1]
    key = chrom + "_" + pos
    type = lline[2]
    if type != "SNV":
        continue
    observations = lline[3]
    FFPE_SNV_artifacts[key] = observations


nr_TMB = 0
header = True
prev_pos = ""
prev_chrom = ""
for line in vcf:
    if header:
        if line[:6] == "#CHROM":
            header = False
        continue
    lline = line.strip().split("\t")
    chrom = lline[0]
    pos = lline[1]
    if chrom == prev_chrom and pos == prev_pos:
        continue
    prev_pos = pos
    prev_chrom = chrom
    key = chrom + "_" + pos
    pos_info[key] = []
    ref = lline[3]
    alt = lline[4]
    filter = lline[6]
    INFO = lline[7]
    if INFO[:3] == "AA=":
        continue
    Variant_type = INFO.split("|")[1].split("&")
    db1000G = INFO.split("|")[41]
    if db1000G == "":
        db1000G = 0
    else:
        db1000G = float(db1000G)
    GnomAD = INFO.split("|")[47]
    if GnomAD == "":
        GnomAD = 0
    else:
        GnomAD = float(GnomAD)
    db = INFO.split("|")[17]
    FORMAT = lline[8].split(":")
    AD_index = 0
    DP_index = 0
    VD_index = 0
    i = 0
    for f in FORMAT:
        if f == "AD":
            AD_index = i
        elif f == "DP":
            DP_index = i
        elif f == "VD":
            VD_index = i
        i += 1
    DATA = lline[9].split(":")
    AD = DATA[AD_index].split(",")
    DP = int(DATA[DP_index])
    VD = int(DATA[VD_index])
    INFO_list = INFO.split(";")
    AF_index = 0
    Caller_index = 0
    i = 0
    for info in INFO_list:
        if info[:3] == "AF=":
            AF_index = i
        if info[:8] == "CALLERS=":
            Caller_index = i
        i += 1
    AF = float(INFO_list[AF_index][3:])
    Callers = INFO_list[Caller_index]
    if Callers.find("vardict") == -1:
        continue

    # Artifact observations
    Observations = 0
    if len(ref) == 1 and len(alt) == 1:
        if key in FFPE_SNV_artifacts:
            Observations = FFPE_SNV_artifacts[key]
    else:
        continue

    # TMB
    if (filter.find("PASS") != -1 and DP > 200 and VD > 20 and AF >= 0.05 and AF <= 0.35 and
            GnomAD <= 0.0001 and db1000G <= 0.0001 and Observations <= 1 and INFO.find("MUC6") == -1 and
            INFO.find("Complex") == -1):
        if ("missense_variant" in Variant_type or
                "splice_region_variant" in Variant_type or
                "splice_acceptor_variant" in Variant_type or
                "stop_gained" in Variant_type or
                # "frameshift_variant" in Variant_type or
                # "protein_altering_variant" in Variant_type or
                "splice_donor_variant" in Variant_type or
                "stop_lost" in Variant_type):
            if len(ref) == 1 and len(alt) == 1:
                nr_TMB += 1

TMB = nr_TMB * 0.78
output_tmb.write("TMB:\t" + str(TMB) + "\n")
output_tmb.write("Variants:\t" + str(nr_TMB) + "\n")
