
#snakemake -p -j 1 --drmaa "-A wp1 -p core -n 1 -t 2:00:00 "  -s ./src/Snakemake/rules/Twist_DNA_yaml/Twist_exome_WP3_yaml.smk

localrules:
    all,
    Create_Twist_DNA_yaml,


rule all:
    input:
        Twist_exome_WP3_yaml="Twist_exome_WP3.yaml",


rule Create_Twist_DNA_yaml:
    input:
        config="Config/Pipeline/config_WP3.yaml",
    output:
        Twist_exome_WP3_yaml="Twist_exome_WP3.yaml",
    run:
        import glob
        import os
        import subprocess

        subprocess.call("cp " + input.config + " " + output.Twist_exome_WP3_yaml, shell=True)
        state = 0
        DNA_sample_list = []
        i = 1
        sample_sheet_name = glob.glob("SampleSheet.csv")
        sample_sheet_name = sample_sheet_name[0]
        infile = open(sample_sheet_name)
        for line in infile:
            if state == 0:
                if line.find("Experiment Name") != -1:
                    state = 1
            elif state == 1:
                if line.find("Sample_ID,Sample_Name") != -1:
                    state = 2
            elif state == 2:
                lline = line.strip().split(",")
                sample = lline[0]
                if sample.find(" ") != -1:
                    print("incorrect sample name: " + sample)
                    quit()
                DNA_sample_list.append([sample, i])
                i += 1
        outfile = open(output.Twist_exome_WP3_yaml, "a")
        outfile.write("DNA_Samples:\n")
        for sample in DNA_sample_list:
            outfile.write("  \"" + sample[0] + "\"\n")
        outfile.close()
