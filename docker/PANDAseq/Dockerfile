FROM ubuntu:16.04
MAINTAINER Inês Mendes <cimendes@medicina.ulisboa.pt>

RUN apt-get update && apt-get -y install build-essential \
    libtool automake zlib1g-dev libbz2-dev pkg-config git

WORKDIR /NGStools

RUN git clone https://github.com/neufeld/pandaseq.git
RUN cd pandaseq && ./autogen.sh && ./configure && make && make install

ENV PATH="/NGStools/pandaseq:$PATH"

