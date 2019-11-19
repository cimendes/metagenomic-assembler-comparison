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
# change this parameter so it will allow reads longer than 128
RUN sed -i 's#static const uint32_t kMaxShortSequence = 128#static const uint32_t kMaxShortSequence = 300#g' /NGStools/idba/src/sequence/short_sequence.h
RUN cd idba && sh build.sh

ENV PATH="/NGStools/idba/bin:$PATH"

WORKDIR /data