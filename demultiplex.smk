# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import src.lib.python.utils as utils

configfile: "demultiplexconfig.yaml"

samples = utils.generate_sample_list_from_samplesheet(config['samplesheet'])

if config.get('lanes', None) is not None:
    fastq1_files = [ "fastq/" + samples[sample]['projectpath'] + sample + "_S" + str(samples[sample]['counter']) + "_" + lane + "_R1_001.fastq.gz" for sample in samples for lane in config['lanes']]
    fastq2_files = [ "fastq/" + samples[sample]['projectpath'] + sample + "_S" + str(samples[sample]['counter']) + "_" + lane + "_R2_001.fastq.gz" for sample in samples for lane in config['lanes']]
else:
    fastq1_files = [ "fastq/" + samples[sample]['projectpath'] + sample + "_S" + str(samples[sample]['counter']) + "_R1_001.fastq.gz" for sample in samples]
    fastq2_files = [ "fastq/" + samples[sample]['projectpath'] + sample + "_S" + str(samples[sample]['counter']) + "_R2_001.fastq.gz" for sample in samples]

rule all:
    input:
        fastq1_files + fastq2_files


include: "src/Snakemake/rules/Fastq/demultiplex.smk"

onsuccess:
    print("Demultiplex workflow finished, no error")
    if config["notification_mail"]:
        shell("mail -s \"Demultiplex workflow completed! \" " + config["notification_mail"] + " < {log}")

onerror:
    print("An error occurred")
    if config["notification_mail"]:
        shell("mail -s \"an error occurred\" " + config["notification_mail"] + " < {log}")
