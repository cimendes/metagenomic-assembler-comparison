FROM ubuntu:16.04
MAINTAINER Inês Mendes <cimendes@medicina.ulisboa.pt>

RUN  apt-get update && apt-get -y install \
    wget \
    lbzip2

WORKDIR /NGStools

RUN wget https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 && \
    tar xf minimap2-2.17_x64-linux.tar.bz2 && rm minimap2-2.17_x64-linux.tar.bz2

ENV PATH="/NGStools/minimap2-2.17_x64-linux:$PATH"
