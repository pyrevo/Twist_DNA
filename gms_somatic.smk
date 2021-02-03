# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import pandas as pd

configfile: "config.yaml"

samples = pd.read_table(config["samples"], index_col="sample")

wildcard_constraints:
    unit="[A-Za-z0-9]+",
    sample="[^.]+",

rule all:
    input:
<<<<<<< HEAD
        [caller + "/" + sample.Index + "." + caller + type for sample in samples.itertuples() for type in [".okAF.vcf"] for caller in ['mutect2', 'vardict', "freebayes"]]
=======
        [caller + "/" + sample.Index + "." + caller  + type for sample in samples.itertuples() for type in [".okAF.vcf.gz"] for caller in ["mutect2", "freebayes"]]
>>>>>>> 39766fb... Some final tweaks.


include: "src/Snakemake/workflow/gms_somatic_workflow.smk"

onsuccess:
    print("Workflow finished, no error")
    if config["notification_mail"]:
        shell("mail -s \"Basic workflow completed! \" " + config["notification_mail"] + " < {log}")

onerror:
    print("An error occurred")
    if config["notification_mail"]:
        shell("mail -s \"an error occurred\" " + config["notification_mail"] + " < {log}")