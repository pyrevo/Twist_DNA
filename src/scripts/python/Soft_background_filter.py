
in_vcf = open(snakemake.input.vcf)
background_panel_filename = snakemake.params.background_panel
out_vcf = open(snakemake.output.vcf, "w")


background_panel_dict = {}
if background_panel_filename != "" :
    background_panel = open(background_panel_filename)
    next(background_panel)
    for line in background_panel:
        lline = line.strip().split("\t")
        chrom = "chr" + lline[0]
        pos = lline[1]
        median = float(lline[2])
        sd = float(lline[3])
        background_panel_dict[chrom + "_" + pos] = [median, sd]


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
    INFO = lline[7]
    if INFO[:3] == "AA=":
        continue
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
    nr_SD = 1000
    if len(ref) == 1 and len(alt) == 1:
        if key in background_panel_dict:
            if background_panel_dict[key][0] > 0.0 :
                nr_SD = (AF - background_panel_dict[key][0]) / background_panel_dict[key][0]
    if nr_SD < 10.0 :
        if filter == "PASS" :
            filter = "LownrSD=" + str(nr_SD)
        else :
            filter += "|LownrSD=" + str(nr_SD)
        lline[6] = filter
        out_vcf.write(lline[0])
        for l in lline[1:] :
            out_vcf.write("\t" + l)
        out_vcf.write("\n")
    else :
        out_vcf.write(line)
out_vcf.close()
