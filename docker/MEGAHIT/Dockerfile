FROM ubuntu:16.04
MAINTAINER C I Mendes, cimendes@medicina.ulisboa.pt

RUN apt-get update

WORKDIR /NGStools

# -- General Dependencies ---
RUN apt-get install -y wget python python3 zlib1g-dev build-essential

# -- megahit ---
RUN wget https://github.com/voutcn/megahit/releases/download/v1.2.8/MEGAHIT-1.2.8-Linux-x86_64-static.tar.gz
RUN tar -xzf MEGAHIT-1.2.8-Linux-x86_64-static.tar.gz && rm MEGAHIT-1.2.8-Linux-x86_64-static.tar.gz
RUN mv MEGAHIT-1.2.8-Linux-x86_64-static megahit

ENV PATH="/NGStools/megahit:$PATH"
