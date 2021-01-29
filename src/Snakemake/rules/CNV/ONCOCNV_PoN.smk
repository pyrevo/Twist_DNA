
configfile: "Twist_DNA.yaml"


rule all:
    input:
        PoN1="DATA/ONCOCNV_Twist_PoN.txt",
        PoN2="DATA/ONCOCNV_Twist_PoN.Processed.txt",


rule fix_bed_file:
    input:
        bed=config["bed"]["bedfile"],
    output:
        bed="CNV/bed/ONCOCNV.bed",
    log:
        "logs/CNV_ONCOCNV/fix_bed_file.log",
    shell:
        "awk 'BEGIN{{ OFS=\"\t\"}}{{ print $1, $2, $3, NR, \"0\", $4 }}' {input.bed} > {output.bed}"


rule Target_GC:
    input:
        bed="CNV/bed/target.bed",
        ref=config["reference"]["ref"],
    output:
        GC="CNV/ONCOCNV_stats/target.GC.txt",
    singularity:
        config["singularity"]["ONCOCNV"]
    shell:
        "perl ONCOCNV/createTargetGC.pl "
        "-bed {input.bed} "
        "-fi {input.ref} "
        "-od CNV/ONCOCNV_stats/ "
        "-of {output.GC}"


rule Normal_levels:
    input:
        #bams=expand("{normal_sample}", normal_sample=config["Normal_samples"]),
        bams=["DNA_bam/" + s + "-ready.bam" for s in config["DNA_Samples"]],
        bed="CNV/bed/ONCOCNV.bed",
    output:
        stats="DATA/ONCOCNV_Twist_PoN.txt",
    singularity:
        config["singularity"]["ONCOCNV"]
    shell:
        # run:
        #     import os
        #     command_line="singularity exec -B /projects/ -B /gluster-storage-volume/ /projects/wp4/nobackup/workspace/jonas_test/ONCOCNV.simg "
        #     command_line += "perl ONCOCNV/ONCOCNV_getCounts.pl getControlStats "
        #     command_line += "-m Ampli "
        #     command_line += "-b " + input[-1]
        #     command_line += " -c \"" + ",".join(input[:-1]) + "\" "
        #     command_line += "-o " + output[0]
        #     print(command_line)
        #     os.system(command_line)
        "perl ONCOCNV/ONCOCNV_getCounts.pl getControlStats -m Ampli -b {input.bed} "
        "-c \"" + ",".join(input[:-1]) + "\" -o {output.stats}"


rule Controls_calls:
    input:
        stats="DATA/ONCOCNV_Twist_PoN.txt",
        GC="CNV/ONCOCNV_stats/target.GC.txt",
    output:
        stats="DATA/CONCOCNV_Twist_PoN.Processed.txt",
    shell:
        "cat ONCOCNV/processControl.R | R --slave "
        "--args {input.stats} {output.stats} {input.GC}"
