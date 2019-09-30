FROM ubuntu:16.04
MAINTAINER Inês Mendes <cimendes@medicina.ulisboa.pt>

RUN apt-get update && apt-get -y install build-essential \
    wget \
    zlib1g-dev \
    perl \
    zip

WORKDIR NGStools/

# nucmer
RUN wget https://github.com/mummer4/mummer/releases/download/v4.0.0beta2/mummer-4.0.0beta2.tar.gz
RUN tar xf mummer-4.0.0beta2.tar.gz && rm mummer-4.0.0beta2.tar.gz
RUN cd mummer-4.0.0beta2 && ./configure && make && make install

# bowtie2
RUN wget https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.3.2/bowtie2-2.3.2-linux-x86_64.zip
RUN unzip bowtie2-2.3.2-linux-x86_64.zip
RUN rm -rf bowtie2-2.3.2-linux-x86_64.zip
ENV PATH="/NGStools/bowtie2-2.3.2:$PATH"

# Samtools
RUN wget https://sourceforge.net/projects/samtools/files/samtools/1.4.1/samtools-1.4.1.tar.bz2
RUN mkdir samtools && tar jxf samtools-1.4.1.tar.bz2 -C samtools --strip-components=1
RUN cd samtools && ./configure --without-curses --disable-bz2 --disable-lzma
RUN cd samtools && make
ENV PATH="/NGStools/samtools:$PATH"


RUN wget https://sourceforge.net/projects/metassembler/files/v1.5/Metassembler.1.5.tar.gz
RUN tar xf Metassembler.1.5.tar.gz && rm Metassembler.1.5.tar.gz
RUN cd Metassembler && make