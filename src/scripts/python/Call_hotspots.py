
input_hotspots = open("DATA/Hotspots_combined.csv")
input_gVCF = open("/projects/wp1/nobackup/ngs/utveckling/Twist_DNA_DATA/gVCF/noUMI/PVAL-35.mutect2.gvcf")

header = True
call_positions = []
for line in input_hotspots :
    if header :
        header = False
        continue
    lline = line.strip().split("\t")
    chrom = str(int(lline[0].split("_")[1].split(".")[0]))
    start_pos = lline[1]
    #end_pos = int(lline[2])
    #gene = lline[3]
    type = lline[6]
    if type != "hotspot" :
        continue
    call_positions.append([chrom, start_pos])
