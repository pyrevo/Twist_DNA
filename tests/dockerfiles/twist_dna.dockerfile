FROM python:3.8

RUN mkdir /Twist_DNA

WORKDIR /Twist_DNA

ADD ./ .

RUN cp -r tests/data /data

RUN cp -r tests/workflow_dry_run/gms_somatic /data_gms_somatic

RUN cp -r tests/workflow_dry_run/demultiplex /data_demultiplex


RUN cp -r tests/workflow_dry_run/twist_dna /data_twist_dna_fgbio
RUN cp tests/workflow_dry_run/twist_dna/configs/Twist_DNA_dup_fgbio.yaml /data_twist_dna_fgbio/Twist_DNA.yaml

RUN cp -r tests/workflow_dry_run/twist_dna /data_twist_dna_markdup
RUN cp tests/workflow_dry_run/twist_dna/configs/Twist_DNA_dup_markduplicates.yaml /data_twist_dna_markdup/Twist_DNA.yaml

RUN cp -r tests/workflow_dry_run/twist_dna /data_twist_dna_gpu
RUN cp tests/workflow_dry_run/twist_dna/configs/Twist_DNA_gpu.yaml  /data_twist_dna_gpu/Twist_DNA.yaml

RUN cp -r tests/workflow_dry_run/twist_dna /data_twist_dna_cutadapt
RUN cp tests/workflow_dry_run/twist_dna/configs/Twist_DNA_cutadapt.yaml  /data_twist_dna_cutadapt/Twist_DNA.yaml

RUN cp -r tests/workflow_dry_run/twist_dna /data_twist_dna_fastp
RUN cp tests/workflow_dry_run/twist_dna/configs/Twist_DNA_fastp.yaml  /data_twist_dna_fastp/Twist_DNA.yaml



RUN pip install -r requirements.txt
