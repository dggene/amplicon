FROM ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive

COPY docker/sources.list /etc/apt/sources.list

# install base require
RUN apt-get update && apt-get install -y --no-install-recommends build-essential \
    python3.6 python3-pip python3-setuptools python3-dev \
    openjdk-8-jdk \
    wget zlib1g.dev curl \
    procps libcurl4-openssl-dev libssl-dev libxml2-dev gnupg2 libbz2-dev liblzma-dev && \
    apt clean && \
    rm -rf /var/lib/apt/lists/*

# install r
COPY docker/r-tuna.list /etc/apt/sources.list.d/r-tuna.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN apt-get update && apt-get install -y --no-install-recommends \
    r-base  && \
    apt clean && \
    rm -rf /var/lib/apt/lists/*

# install git
RUN apt-get update && apt-get install -y \
    git  && \
    apt clean && \
    rm -rf /var/lib/apt/lists/*

# begin BWA
ENV BWA_VERSION 0.7.5a
RUN wget https://github.com/lh3/bwa/archive/${BWA_VERSION}.tar.gz -O /tmp/${BWA_VERSION}.tar.gz && \
	tar -C /tmp -zxvf /tmp/${BWA_VERSION}.tar.gz  && \
    cd /tmp/bwa-${BWA_VERSION} && \
    make clean all && \
    cp -rf bwa /usr/local/bin/ && \
    rm -rf /tmp/${BWA_VERSION}.tar.gz /tmp/bwa-${BWA_VERSION}
# end BWA

# begin soapnuke
ENV SOAPNUKE_VERSION 2.1.0
RUN wget https://github.com/BGI-flexlab/SOAPnuke/archive/SOAPnuke${SOAPNUKE_VERSION}.tar.gz -O /tmp/SOAPnuke${SOAPNUKE_VERSION}.tar.gz && \
	tar -C /tmp -zxvf /tmp/SOAPnuke${SOAPNUKE_VERSION}.tar.gz  && \
    cd /tmp/SOAPnuke-SOAPnuke${SOAPNUKE_VERSION} && \
    make && \
    cp -rf SOAPnuke /usr/local/bin/ && \
    rm -rf /tmp/SOAPnuke${SOAPNUKE_VERSION}.tar.gz /tmp/SOAPnuke-SOAPnuke${SOAPNUKE_VERSION} 
#end soapnuke

RUN mkdir -p /usr/picard /usr/gatk
RUN wget https://gcs.obs.cn-north-4.myhuaweicloud.com/tools/picard.jar  -O /usr/picard/picard.jar

ENV GATK_VERSION 3.7-0-gcfedb67
RUN wget https://gcs.obs.cn-north-4.myhuaweicloud.com/tools/GenomeAnalysisTK-${GATK_VERSION}.tar.bz2 -O /tmp/GenomeAnalysisTK-${GATK_VERSION}.tar.bz2 && \
    tar -C /tmp -xjf /tmp/GenomeAnalysisTK-${GATK_VERSION}.tar.bz2 && \
    cp -rf /tmp/GenomeAnalysisTK-${GATK_VERSION}/GenomeAnalysisTK.jar /usr/gatk/gatk.jar && \
    rm -rf /tmp/GenomeAnalysisTK-${GATK_VERSION}.tar.bz2 /tmp/GenomeAnalysisTK-${GATK_VERSION}

ADD docker/package.r /tmp/package.r
RUN Rscript /tmp/package.r

COPY docker/requirements.txt /tmp/requirements.txt
RUN pip3 install -r /tmp/requirements.txt  --index-url=https://mirrors.aliyun.com/pypi/simple/ --trusted-host=mirrors.aliyun.com

ENV NXF_OFFLINE="TRUE"
ENV NEXTFLOW_VERSION 20.10.0
# https://github.com/nextflow-io/nextflow/releases/download/v${NEXTFLOW_VERSION}/nextflow-${NEXTFLOW_VERSION}-all
RUN wget https://gcs.obs.cn-north-4.myhuaweicloud.com:443/tools/nextflow-${NEXTFLOW_VERSION}-all -O nextflow && \
    mv nextflow /usr/local/bin && \
    chmod 755 /usr/local/bin/nextflow

WORKDIR /code

ADD . /code/

WORKDIR /home/work_dir
CMD /bin/bash