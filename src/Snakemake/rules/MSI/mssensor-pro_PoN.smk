rule all:
    input:
        PoN_list="/projects/wp1/nobackup/ngs/utveckling/Twist_DNA_DATA/MSI_PoN_new/Msisensor_pro_reference.list_baseline",


rule msisensor_pro_scan:
    input:
        ref=config["reference"]["ref"],
    output:
        PoN_list="/projects/wp1/nobackup/ngs/utveckling/Twist_DNA_DATA/MSI_PoN_new/Msisensor_pro_reference.list",
    container:
        config["singularity"].get("msisensor-pro", config["singularity"].get("default", ""))
    shell:
        "msisensor-pro scan -d {input.ref} -o {output.PoN_list}"


rule msisensor_pro_baseline:
    input:
        Normal_bams="/projects/wp1/nobackup/ngs/utveckling/Twist_DNA_DATA/MSI_PoN_new/configure.txt",
        PoN=config["PoN"]["msisensor-pro"],
    output:
        PoN_list="/projects/wp1/nobackup/ngs/utveckling/Twist_DNA_DATA/MSI_PoN_new/Msisensor_pro_reference.list_baseline",
    params:
        extra="-c 50",  # -c = minimal coverage
        out_dir="/projects/wp1/nobackup/ngs/utveckling/Twist_DNA_DATA/MSI_PoN_new/",
        out_file="/projects/wp1/nobackup/ngs/utveckling/Twist_DNA_DATA/MSI_PoN_new/reference.list_baseline",
    container:
        config["singularity"].get("msisensor-pro", config["singularity"].get("default", ""))
    shell:
        "msisensor-pro baseline {params.extra} -d {output.PoN_list} -i {input.Normal_bams} -o {params.out_dir} && "
        "mv {params.out_file} {output.PoN_list}"
