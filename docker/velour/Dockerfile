FROM ubuntu:16.04
MAINTAINER Inês Mendes <cimendes@medicina.ulisboa.pt>

RUN apt-get update && apt-get -y install build-essential \
    git \
    zlib1g-dev \
    wget

WORKDIR /NGStools

RUN git clone https://github.com/jjcook/velour.git
RUN cd velour && ls velour_minikmer_ptables*.tar.bz2 | xargs -n 1 tar xjvf && make
ENV PATH="/NGStools/velour:$PATH"
