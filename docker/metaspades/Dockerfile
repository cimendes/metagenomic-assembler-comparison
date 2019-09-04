FROM ubuntu:16.04
MAINTAINER Inês Mendes <cimendes@medicina.ulisboa.pt>

RUN  apt-get update && apt-get -y install git \
    wget \
    python \
    build-essential \
    cmake

RUN apt-get -y install zlib1g-dev libbz2-dev

WORKDIR /NGStools

RUN wget https://github.com/ablab/spades/archive/spaligner-paper.tar.gz
RUN tar xf spaligner-paper.tar.gz && rm spaligner-paper.tar.gz
RUN mv spades-spaligner-paper spades
RUN cd spades/assembler && sh spades_compile.sh
ENV PATH="/NGStools/spades/assembler:$PATH"
