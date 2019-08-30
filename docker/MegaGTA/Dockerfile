FROM ubuntu:16.04
MAINTAINER Inês Mendes <cimendes@medicina.ulisboa.pt>

RUN  apt-get update && apt-get -y install git \
    g++ \
    wget \
    make \
    ant \
    zlib1g
RUN apt-get -y install openjdk-8-jdk

WORKDIR /NGStools

# HMMER 3.1
RUN wget http://eddylab.org/software/hmmer/hmmer-3.1b2.tar.gz && tar xf hmmer-3.1b2.tar.gz  && rm hmmer-3.1b2.tar.gz
RUN cd hmmer-3.1b2/ && sh configure && make && make check && make install
RUN cd hmmer-3.1b2/easel; make install

# UCHIME
RUN wget http://drive5.com/uchime/uchime4.2.40_src.tar.gz && tar xf uchime4.2.40_src.tar.gz && rm uchime4.2.40_src.tar.gz
RUN cd uchime4.2.40_src && make
ENV PATH="/NGStools/uchime4.2.40_src:$PATH"

# MegaGTA
RUN git clone https://github.com/HKU-BAL/megagta.git
RUN cd megagta && git submodule update --init --recursive
RUN cd megagta && ./make.sh
RUN sed -i 's#UCHIME=/nas5/ykhuang/uchime4.2.40_i86linux32 #UCHIME=/NGStools/uchime4.2.40_src #g' /NGStools/megagta/bin/post_proc.sh
RUN sed -i 's#HMMALIGN=/nas5/ykhuang/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmalign #HMMALIGN=/NGStools/uchime4.2.40_src #g' /NGStools/megagta/bin/post_proc.sh

