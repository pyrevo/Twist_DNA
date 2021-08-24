
in_vcf = open(snakemake.input.vcf)
artifacts = open(snakemake.params.artifacts)
out_vcf = open(snakemake.output.vcf, "w")


artifact_dict = {}
next(artifacts)
for line in artifacts:
    lline = line.strip().split("\t")
    chrom = lline[0]
    pos = lline[1]
    type = lline[2]
    observations = int(lline[3])
    artifact_dict[chrom + "_" + pos] = [type, observations]


header = True
for line in in_vcf :
    if header:
        out_vcf.write(line)
        if line[:6] == "#CHROM":
            header = False
        continue
    lline = line.strip().split("\t")
    chrom = lline[0]
    pos = lline[1]
    key = chrom + "_" + pos
    ref = lline[3]
    alt = lline[4]
    filter = lline[6]
    Observations = 0
    if len(ref) == 1 and len(alt) == 1:
        if key in artifact_dict:
            if artifact_dict[key][0] == "SNV" :
                Observations = artifact_dict[key][1]
    else:
        if key in artifact_dict:
            if artifact_dict[key][0] == "INDEL" :
                Observations = artifact_dict[key][1]
    if Observations >= 2 :
        if filter == "PASS" :
            filter = "Arti=" + str(Observations)
        else :
            filter += "|Arti=" + str(Observations)
        lline[6] = filter
        out_vcf.write(lline[0])
        for l in lline[1:] :
            out_vcf.write("\t" + l)
        out_vcf.write("\n")
    else :
        out_vcf.write(line)
out_vcf.close()
