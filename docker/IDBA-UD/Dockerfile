FROM ubuntu:16.04
MAINTAINER Inês Mendes <cimendes@medicina.ulisboa.pt>

RUN  apt-get update && apt-get -y install git \
    gcc \
    autoconf \
    automake \
    g++ \
    make


WORKDIR /NGStools

RUN git clone https://github.com/loneknightpy/idba.git
RUN cd idba && sh build.sh

ENV PATH="/NGStools/idba/bin:$PATH"

WORKDIR /data