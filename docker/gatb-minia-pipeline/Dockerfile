FROM ubuntu:16.04
MAINTAINER Inês Mendes <cimendes@medicina.ulisboa.pt>

RUN apt-get update && apt-get -y install build-essential \
    git \
    python-dev \
    python-pip \
    zlib1g-dev \
    python-scipy

WORKDIR /NGStools

# BWA dependency
RUN git clone https://github.com/lh3/bwa.git
RUN cd bwa && make
ENV PATH="/NGStools/bwa:$PATH"

RUN pip install mathstats networkx==1.9 pyfasta pysam==0.8.3 matplotlib==2.2.4

RUN git clone --recursive https://github.com/GATB/gatb-minia-pipeline
ENV PATH="/NGStools/gatb-minia-pipeline:$PATH"