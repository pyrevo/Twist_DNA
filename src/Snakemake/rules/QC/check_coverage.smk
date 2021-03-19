

rule check_DNA_coverage:
    input:
        bam="DNA_bam/{sample}-ready.bam",
        bai="DNA_bam/{sample}-ready.bam.bai",
        bed="DATA/Hotspots_combined.csv",
        # bediles = "Mutations_Colon_20171219.csv", "Mutations_Lung_20190211.csv", "Mutations_Melanom_20190211.csv", "Mutations_Gist_20190211.csv","Mutations_Ovarial_20170125_header_ok.csv"]
    output:
        coverage="qc/{sample}/Low_coverage_positions.txt",
        coverage2="qc/{sample}/All_coverage_positions.txt",
    log:
        "logs/qc/{sample}_check_DNA_coverage.log",
    container:
        config["singularity"].get("python", config["singularity"].get("default", ""))
    script:
        "../../../scripts/python/check_coverage.py"
