FROM ubuntu:16.04
MAINTAINER Inês Mendes <cimendes@medicina.ulisboa.pt>

RUN apt-get update && apt-get -y install build-essential \
    wget \
    git

WORKDIR /NGStools

# Dependencies
# BCALM
RUN wget https://github.com/GATB/bcalm/releases/download/v2.2.3/bcalm-binaries-v2.2.3-Linux.tar.gz
RUN tar xf bcalm-binaries-v2.2.3-Linux.tar.gz && rm bcalm-binaries-v2.2.3-Linux.tar.gz

ENV PATH="/NGStools/bcalm-binaries-v2.2.3-Linux.tar.gz/bin:${PATH}"

WORKDIR /data
