FROM ubuntu:16.04
MAINTAINER Inês Mendes <cimendes@medicina.ulisboa.pt>

RUN apt-get update && apt-get -y install build-essential \
    wget

WORKDIR /NGStools

# Install dependencies
#    Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN chmod +x Miniconda3-latest-Linux-x86_64.sh
RUN ./Miniconda3-latest-Linux-x86_64.sh -b -p /NGStools/miniconda
ENV PATH="/NGStools/miniconda/bin:$PATH"
#    Configure channels for bioconda
RUN conda config --add channels defaults && conda config --add channels conda-forge && conda config --add channels bioconda

# Install MetaVelvet-SL
RUN conda install -c bioconda metavelvet-sl

# Download pre-trained classification models
RUN http://metavelvet.dna.bio.keio.ac.jp/srcSL/LibraryPretrainedClassification.tar.gz
RUN tar xf LibraryPretrainedClassification.tar.gz && rm LibraryPretrainedClassification.tar.gz

WORKDIR /data