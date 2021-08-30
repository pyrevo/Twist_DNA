
from pysam import VariantFile

in_vcf_filename = snakemake.input.vcf
artifacts = open(snakemake.params.artifacts)
out_vcf_filename = snakemake.output.vcf



artifact_dict = {}
header = True
caller_list = []
for line in artifacts:
    lline = line.strip().split("\t")
    if header:
        caller_list = lline[3:]
        header = False
        continue
    chrom = lline[0]
    pos = lline[1]
    type = lline[2]
    observations = lline[3:]
    artifact_dict[chrom + "_" + pos] = [type, observations]


in_vcf = VariantFile(in_vcf_filename)
new_header = in_vcf.header
new_header.filters.add("Artifact", None, None, "SNV or INDEL observed in other samples")
new_header.info.add(
    "Artifact", "1", "String", "Number of observations of SNV or INDEL in panel samples per caller: " + "|".join(caller_list)
)
out_vcf = VariantFile(out_vcf_filename, 'w', header=new_header)
out_vcf.close()
in_vcf.close()


out_vcf = open(out_vcf_filename, "a")
in_vcf = open(in_vcf_filename)
header = True
for line in in_vcf:
    if header:
        # out_vcf.write(line)
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
    Observations = []
    if len(ref) == 1 and len(alt) == 1:
        if key in artifact_dict:
            if artifact_dict[key][0] == "SNV":
                Observations = artifact_dict[key][1]
    else:
        if key in artifact_dict:
            if artifact_dict[key][0] == "INDEL":
                Observations = artifact_dict[key][1]
    max_observations = 0
    for obs in Observations:
        if int(obs) > max_observations :
            max_observations = int(obs)
    if max_observations >= 2:
        if filter == "PASS":
            filter = "Artifact"
        else:
            filter += ";Artifact"
        lline[6] = filter
    INFO = lline[7]
    Artifact_string = "Artifact=" + Observations[0]
    if len(Observations) > 1:
        for obs in Observations[1:]:
            Artifact_string += "|" + obs
    INFO = Artifact_string + ";" + INFO
    lline[7] = INFO
    out_vcf.write(lline[0])
    for column in lline[1:]:
        out_vcf.write("\t" + column)
    out_vcf.write("\n")
out_vcf.close()
