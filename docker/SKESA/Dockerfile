FROM ubuntu:16.04
MAINTAINER Inês Mendes <cimendes@medicina.ulisboa.pt>

# not working

RUN apt-get update && apt-get -y install \
	bzip2 \
	git \
	python3 \
	wget

WORKDIR /NGStools

RUN wget https://github.com/ncbi/SKESA/releases/download/v2.3.0/skesa.centos6.10
RUN mv skesa.centos6.10 skesa && chmod +x skesa

ENV PATH="/NGStools:$PATH"