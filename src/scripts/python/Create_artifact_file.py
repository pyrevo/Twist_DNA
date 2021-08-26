
in_vcfs = open(snakemake.input.vcfs)
caller_list = snakemake.params.callers
out_artifacts = open(snakemake.output.artifacts, "w")


def FFPE_artifacts():
    ffpe_snv_dict = {}
    ffpe_indel_dict = {}
    for vcf_filename in in_vcfs:
        vcf = open(vcf_filename.strip())
        header = True
        prev_pos = 0
        for line in vcf:
            if header:
                if line[:6] == "#CHROM":
                    header = False
                continue
            lline = line.strip().split("\t")
            chrom = lline[0]
            pos = int(lline[1])
            if prev_pos == pos:
                continue
            key = chrom + "_" + str(pos)
            ref = lline[3]
            alt = lline[4]
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
            prev_pos = pos
            AF = float(INFO_list[AF_index][3:])
            Callers = INFO_list[Caller_index].split("=")[1].split(",")
            if len(ref) == 1 and len(alt) == 1:
                if key not in ffpe_snv_dict:
                    ffpe_snv_dict[key] = {}
                    for caller in caller_list :
                        ffpe_snv_dict[key][caller] = 0
                for caller in Callers :
                    ffpe_snv_dict[key][caller] += 1
            else:
                if key not in ffpe_indel_dict:
                    ffpe_indel_dict[key] = {}
                    for caller in caller_list :
                        ffpe_indel_dict[key][caller] = 0
                for caller in Callers :
                    ffpe_indel_dict[key][caller] += 1
        vcf.close()

    return ffpe_snv_dict, ffpe_indel_dict


def write_artifacts():
    out_artifacts.write("Chrom\tPos\tVariant_type")
    for caller in caller_list :
        out_artifacts.write("\t" + caller)
    out_artifacts.write("\n")
    for key in ffpe_snv_dict:
        out_artifacts.write(key.split("_")[0] + "\t" + key.split("_")[1] + "\tSNV")
        for caller in caller_list :
            out_artifacts.write("\t" + str(ffpe_snv_dict[key][caller]))
        out_artifacts.write("\n")
    for key in ffpe_indel_dict:
        out_artifacts.write(key.split("_")[0] + "\t" + key.split("_")[1] + "\tINDEL")
        for caller in caller_list :
            out_artifacts.write("\t" + str(ffpe_indel_dict[key][caller]))
        out_artifacts.write("\n")


ffpe_snv_dict, ffpe_indel_dict = FFPE_artifacts()
write_artifacts()
