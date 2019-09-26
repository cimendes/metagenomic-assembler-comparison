FROM ubuntu:16.04
MAINTAINER Inês Mendes <cimendes@medicina.ulisboa.pt>

RUN apt-get update && apt-get -y install build-essential \
    wget \
    zlib1g-dev

WORKDIR /NGStools

RUN wget https://www.ebi.ac.uk/~zerbino/velvet/velvet_1.2.10.tgz
RUN tar xf velvet_1.2.10.tgz && rm velvet_1.2.10.tgz

RUN cd velvet_1.2.10 && make
ENV PATH="/NGStools/velvet_1.2.10:$PATH"


