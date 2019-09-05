FROM ubuntu:16.04
MAINTAINER Inês Mendes <cimendes@medicina.ulisboa.pt>

RUN  apt-get update && apt-get -y install git \
    wget \
    gcc \
    make

WORKDIR /NGStools

# INSTALL HMMER 3.0
RUN wget http://eddylab.org/software/hmmer/hmmer-3.1b2.tar.gz && tar xf hmmer-3.1b2.tar.gz  && rm hmmer-3.1b2.tar.gz
RUN cd hmmer-3.1b2/ && sh configure && make && make check && make install
RUN cd hmmer-3.1b2/easel; make install

# INSTALL MINICONDA
RUN wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
RUN chmod +x Miniconda2-latest-Linux-x86_64.sh
RUN ./Miniconda2-latest-Linux-x86_64.sh -b -p /NGStools/miniconda
ENV PATH="/NGStools/miniconda/bin:$PATH"
# Configure channels for bioconda
RUN conda config --add channels defaults && conda config --add channels conda-forge && conda config --add channels bioconda

# INSTALL PYTHON2 and Biopython1.64
RUN conda install python=2.7 biopython

# DOWNLOAD Pfam profile HMMs
RUN wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam30.0/Pfam-A.hmm.gz
RUN gunzip Pfam-A.hmm.gz

# INSTALL Snowball
RUN git clone https://github.com/hzi-bifo/snowball.git
ENV PATH="/NGStools/snowball/algbioi/ga:$PATH"
RUN export PYTHONPATH=/NGStools/snowball:$PYTHONPATH