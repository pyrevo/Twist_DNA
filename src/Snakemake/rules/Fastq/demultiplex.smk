
fastq1_files = []
fastq2_files = []

if config["DNA_Samples"] != "No DNA" :
    S_dna = []
    for s in config["DNA_Samples"].values() :
        S_dna.append(s)
    fastq1_files += ["fastq_temp/DNA/" + s + "_" + i + "_R1_001.fastq.gz" for s,i in zip(config["DNA_Samples"], S_dna)]
    fastq2_files +=  ["fastq_temp/DNA/" + s + "_" + i + "_R2_001.fastq.gz" for s,i in zip(config["DNA_Samples"], S_dna)]



rule demultiplex:
    output:
        fastq1 = fastq1_files,
        fastq2 = fastq2_files
    params:
        runfolder = config["Runfolder"],
        sample_sheet = config["Sample_sheet"]
    run:
        shell("module add bcl2fastq/2.17.1.14; bcl2fastq -i {params.runfolder}/Data/Intensities/BaseCalls -o fastq_temp/ --sample-sheet {params.sample_sheet} --barcode-mismatches 1 --no-lane-splitting -r 16 -d 16 -p 16")
