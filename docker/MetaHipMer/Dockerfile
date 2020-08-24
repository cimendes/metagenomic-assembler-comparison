FROM ubuntu:16.04
MAINTAINER Inês Mendes <cimendes@medicina.ulisboa.pt>

RUN  apt-get update && apt-get -y install \
    wget \
    build-essential \
    gcc

WORKDIR /NGStools

RUN wget https://sourceforge.net/projects/hipmer/files/HipMer-1.2.2.tar.gz && \
    tar xf HipMer-1.2.2.tar.gz && rm HipMer-1.2.2.tar.gz

RUN mv HipMer-* HipMer

WORKDIR HipMer

RUN mkdir build && cd build && cmake .. && make install

ENV PATH="/NGStools/HipMer:$PATH"
