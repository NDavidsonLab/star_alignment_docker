# FROM mgibio/star:2.7.0f

FROM ubuntu:22.04

RUN apt update && \
    apt install -y wget curl build-essential ca-certificates python3-pip unzip

# install fastp
RUN wget http://opengene.org/fastp/fastp.1.3.3 && \
    mv fastp.1.3.3 /usr/local/sbin/fastp && \
    chmod a+x /usr/local/sbin/fastp

# install salmon
RUN cd /tmp && \
    wget https://github.com/COMBINE-lab/salmon/releases/download/v1.11.4/salmon-linux-x86_64.tar.gz && \
    tar zxf salmon-linux-x86_64.tar.gz && \
    mv salmon-linux-x86_64/bin/salmon /usr/local/sbin/salmon && \
    chmod +x /usr/local/sbin/salmon

# install openjdk, a requirement for fastqc
RUN apt install -y --no-install-recommends openjdk-17-jre-headless

# install fastqc via download
RUN cd /tmp && \
    wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && \
    unzip fastqc_v0.12.1.zip && \
    ln -s /tmp/FastQC/fastqc /usr/local/sbin/fastqc && \
    chmod +x /usr/local/sbin/fastqc

ENV CLASSPATH=/tmp/FastQC/fastqc

# install STAR
RUN wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.11b.tar.gz \
    && tar -xzf 2.7.11b.tar.gz \
    && mv STAR-2.7.11b/bin/Linux_x86_64_static/STAR /usr/local/bin/STAR \
    && rm -rf STAR-2.7.11b*

# install multiqc
RUN pip3 install multiqc

# install fastp
RUN apt install -y fastp

# install gffread
RUN apt install -y gffread
