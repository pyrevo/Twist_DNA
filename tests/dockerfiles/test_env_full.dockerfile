FROM jenkins/inbound-agent
ARG DEBIAN_FRONTEND=noninteractive

USER root

RUN wget -O- http://neuro.debian.net/lists/xenial.us-ca.full | tee /etc/apt/sources.list.d/neurodebian.sources.list
RUN apt-key adv --recv-keys --keyserver hkp://pool.sks-keyservers.net:80 0xA5D32F012649A5A9
RUN apt-get -qy update

RUN apt-get install -qy singularity-container python3 python3-dev python3-pip mailutils  gawk vim

ADD requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

RUN mkdir /singularity
WORKDIR /singularity

RUN singularity pull docker://gmsuppsala/somatic:develop
RUN singularity pull docker://broadinstitute/gatk:4.1.9.0
RUN singularity pull docker://rjmashl/msisensor2_mgi:0.1
RUN singularity pull docker://ensemblorg/ensembl-vep:release_99.0
