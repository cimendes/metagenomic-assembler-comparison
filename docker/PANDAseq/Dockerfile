FROM ubuntu:16.04
MAINTAINER Inês Mendes <cimendes@medicina.ulisboa.pt>

RUN apt-get update && apt-get -y install build-essential \
    libtool automake zlib1g-dev libbz2-dev pkg-config git wget

WORKDIR /NGStools

RUN wget https://github.com/neufeld/pandaseq/archive/v2.11.tar.gz
RUN tar xf v2.11.tar.gz && rm v2.11.tar.gz
RUN cd pandaseq-2.11 && ./autogen.sh && ./configure && make && make install
RUN mv pandaseq-2.11 pandaseq
RUN ldconfig

ENV PATH="/NGStools/pandaseq:$PATH"

