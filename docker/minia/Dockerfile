FROM ubuntu:16.04
MAINTAINER Inês Mendes <cimendes@medicina.ulisboa.pt>

RUN apt-get update && apt-get -y install build-essential \
    zlib1g-dev \
    wget \
    python-dev

WORKDIR /usr/local/src

# Install updated cmake
RUN wget https://github.com/Kitware/CMake/releases/download/v3.18.0/cmake-3.18.0.tar.gz
RUN tar xf cmake-3.18.0.tar.gz && rm cmake-3.18.0.tar.gz
RUN cd cmake-3.18.0 && ./bootstrap && make && make install

WORKDIR /NGStools
RUN wget https://github.com/GATB/minia/releases/download/v3.2.4/minia-v3.2.4-Source.tar.gz
RUN tar xf minia-v3.2.4-Source.tar.gz && rm minia-v3.2.4-Source.tar.gz
RUN mv minia-v3.2.4-Source minia
RUN cd minia && sh INSTALL
ENV PATH="/NGStools/minia/build/bin:$PATH"

WORKDIR /data
