FROM ubuntu:16.04
MAINTAINER Inês Mendes <cimendes@medicina.ulisboa.pt>

RUN  apt-get update && apt-get -y install git \
    wget \
    python \
    build-essential \
    cmake

RUN apt-get -y install zlib1g-dev libbz2-dev

WORKDIR /NGStools

RUN wget https://github.com/ablab/spades/releases/download/v3.14.0/SPAdes-3.14.0.tar.gz
RUN tar xf SPAdes-3.14.0.tar.gz && rm SPAdes-3.14.0.tar.gz
RUN mv SPAdes-3.14.0 spades
ENV PATH="/NGStools/spades/bin:$PATH"
