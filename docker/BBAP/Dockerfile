FROM ubuntu:16.04
MAINTAINER Inês Mendes <cimendes@medicina.ulisboa.pt>

RUN apt-get update && apt-get -y install build-essential \
    wget \
    unrar \
    perl

WORKDIR /NGStools

# Dependencies
## BLAST
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-x64-linux.tar.gz && tar xf ncbi-blast-2.9.0+-x64-linux.tar.gz && rm ncbi-blast-2.9.0+-x64-linux.tar.gz
ENV PATH="/NGStools/ncbi-blast-2.9.0+/bin:${PATH}"

# BBAP
RUN wget http://homepage.ntu.edu.tw/~youylin/BBAP.rar
RUN unrar x BBAP.rar && rm BBAP.rar

RUN chmod +x /NGStools/BBAP/*.pl
ENV PATH="/NGStools/BBAP:${PATH}"

WORKDIR /data




