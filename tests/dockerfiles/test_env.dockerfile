FROM centos:8

RUN yum update -y && \
    yum install -y epel-release && \
    yum update -y && \
    yum install -y singularity-runtime singularity mailx  && \
    yum install -y python38-devel python38-pip

ADD requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

RUN mkdir /singularity
WORKDIR /singularity

RUN singularity pull docker://gmsuppsala/somatic:develop
RUN singularity pull docker://broadinstitute/gatk:4.1.9.0
RUN singularity pull docker://rjmashl/msisensor2_mgi:0.1
RUN singularity pull docker://ensemblorg/ensembl-vep:release_99.0
