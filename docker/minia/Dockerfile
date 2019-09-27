FROM ubuntu:16.04
MAINTAINER Inês Mendes <cimendes@medicina.ulisboa.pt>

RUN apt-get update && apt-get -y install build-essential \
    git \
    zlib1g-dev \
    wget \
    python-dev

WORKDIR /usr/local/src

# Install updated cmake
RUN wget https://github.com/Kitware/CMake/releases/download/v3.15.3/cmake-3.15.3.tar.gz
RUN tar xf cmake-3.15.3.tar.gz && rm cmake-3.15.3.tar.gz
RUN cd cmake-3.15.3 && ./bootstrap && make && make install

WORKDIR /NGStools
RUN git clone --recursive https://github.com/GATB/minia.git
RUN cd minia && sh INSTALL
ENV PATH="/NGStools/minia:$PATH"