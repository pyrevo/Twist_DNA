
vcf = open(snakemake.input.vcf)
artifacts = open(snakemake.input.artifacts)
output_tmb = open(snakemake.input.tmb, "w")


ffpe_snv_dict, ffpe_indel_dict = FFPE_artifacts()



nr_TMB = 0
header = True
prev_pos = ""
prev_chrom = ""
for line in vcf :
    if header :
        if line[:6] == "#CHROM" :
            header = False
        continue
    lline = line.strip().split("\t")
    chrom = lline[0]
    pos = lline[1]
    if chrom == prev_chrom and pos == prev_pos :
        continue
    prev_pos = pos
    prev_chrom = chrom
    key = chrom + "_" + pos
    pos_info[key] = []
    ref = lline[3]
    alt = lline[4]
    filter = lline[6]
    INFO = lline[7]
    if INFO[:3] == "AA=" :
        continue
    Variant_type = INFO.split("|")[1].split("&")
    db1000G = INFO.split("|")[41]
    if db1000G == "" :
        db1000G = 0
    else :
        db1000G = float(db1000G)
    GnomAD = INFO.split("|")[47]
    if GnomAD == "" :
        GnomAD = 0
    else :
        GnomAD = float(GnomAD)
    db = INFO.split("|")[17]
    FORMAT = lline[8].split(":")
    AD_index = 0
    DP_index = 0
    VD_index = 0
    i = 0
    for f in FORMAT :
        if f == "AD" :
            AD_index = i
        elif f == "DP" :
            DP_index = i
        elif f == "VD" :
            VD_index = i
        i += 1
    DATA = lline[9].split(":")
    AD = DATA[AD_index].split(",")
    DP = int(DATA[DP_index])
    VD = int(DATA[VD_index])
    INFO_list = INFO.split(";")
    AF_index = 0
    i = 0
    for info in INFO_list:
        if info[:3] == "AF=" :
            AF_index = i
        i += 1
    AF = float(INFO_list[AF_index][3:])

    #SNVs and INDELs and artifacts
    Observations = 0
    if len(ref) == 1 and len(alt) == 1 :
        if key in ffpe_snv_dict :
            Observations = ffpe_snv_dict[key]
    elif len(ref) > 1 :
        if key in ffpe_indel_dict :
            Observations = ffpe_indel_dict[key]
    elif len(alt) > 1 :
        if key in ffpe_indel_dict :
            Observations = ffpe_indel_dict[key]


    #TMB
    if filter.find("PASS") != -1 and DP > 200 and VD > 20 and AF >= 0.05 and AF <= 0.35 and GnomAD <= 0.0001 and db1000G <= 0.0001 and Observations <= 2 and INFO.find("MUC6") == -1 and INFO.find("Complex") == -1 :
        if (#"frameshift_variant" in Variant_type or
            #"inframe_deletion" in Variant_type or
            "missense_variant" in Variant_type or
            "splice_region_variant" in Variant_type or
            "splice_acceptor_variant" in Variant_type or
            "stop_gained" in Variant_type or
            #"protein_altering_variant" in Variant_type or
            "splice_donor_variant" in Variant_type or
            #"inframe_insertion" in Variant_type or
            "stop_lost" in Variant_type) :
            #print(line.strip())
            #print(Variant_type)
            if len(ref) == 1 and len(alt) == 1 :
                nr_TMB += 1

TMB = nr_TMB * 0.75
output_tmb.write("TMB:\t" + str(TMB) + "\n")
output_tmb.write("Variants:\t" + str(nr_TMB) + "\n")
