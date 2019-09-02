FROM ubuntu:16.04
MAINTAINER Inês Mendes <cimendes@medicina.ulisboa.pt>

RUN  apt-get update && apt-get -y install git \
    wget \
    python

WORKDIR /NGStools

RUN wget https://github.com/ablab/spades/archive/spaligner-paper.tar.gz
RUN tar xf spaligner-paper.tar.gz && rm spaligner-paper.tar.gz
RUN mv spades-spaligner-paper spades

ENV PATH="/NGStools/spades/assembler:$PATH"
