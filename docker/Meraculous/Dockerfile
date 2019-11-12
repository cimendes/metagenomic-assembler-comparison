FROM ubuntu:16.04
MAINTAINER Inês Mendes <cimendes@medicina.ulisboa.pt>

RUN apt-get update && apt-get -y install build-essential \
    wget \
    libboost-all-dev \
    perl

RUN whereis boost

WORKDIR /NGStools

RUN wget ftp://ftp.jgi-psf.org/pub/JGI_data/meraculous/Meraculous-1.4.4-rel.tar
RUN tar xf Meraculous-1.4.4-rel.tar && rm Meraculous-1.4.4-rel.tar

RUN sed -i 's#g++ -fPIC -c KeylessHash_wrap.cxx -I /usr/lib/perl/5.10/CORE/ -I ${BOOST_DIR} #g++ -fPIC -c KeylessHash_wrap.cxx -I /usr/lib/perl/5.10/CORE/ -I /usr/include/boost #g' /NGStools/Meraculous-1.4.4-rel/MERACULOUS/build_swig
RUN sed -i 's#g++ -I ${BOOST_DIR} -fPIC -c KeylessHash.cpp #g++ -I /usr/include/boost -fPIC -c KeylessHash.cpp #g' /NGStools/Meraculous-1.4.4-rel/MERACULOUS/build_swig

RUN cd /NGStools/Meraculous-1.4.4-rel/MERACULOUS/ && ./build_swig
