FROM ubuntu:16.04
MAINTAINER Inês Mendes <cimendes@medicina.ulisboa.pt>

RUN apt-get update && apt-get -y install build-essential \
    git \
    zlib1g-dev \
    wget \
    perl \
    bioperl

WORKDIR /NGStools

# Install velvet
RUN wget https://github.com/adigenova/wengan/releases/download/0.1/wengan-v0.1-bin-Linux.tar.gz
RUN tar xf wengan-v0.1-bin-Linux.tar.gz && rm wengan-v0.1-bin-Linux.tar.gz
RUN cd wengan-v0.1-bin-Linux && ls
ENV PATH="/NGStools/wengan-v0.1-bin-Linux:$PATH"
