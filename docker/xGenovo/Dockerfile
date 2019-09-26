FROM ubuntu:16.04
MAINTAINER Inês Mendes <cimendes@medicina.ulisboa.pt>

RUN apt-get update && apt-get -y install build-essential \
    wget \
    python \
    python-pip \
    libgsl2

RUN pip install gdown

RUN ln -s /usr/lib/x86_64-linux-gnu/libgsl.so.19 /usr/lib/x86_64-linux-gnu/libgsl.so.0

WORKDIR /NGStools

RUN gdown https://drive.google.com/uc?id=0B48F8o_nwW-nTTliZ3RlWVFhOEU
RUN tar xf Xgenovo.tar.gz && rm Xgenovo.tar.gz

RUN chmod +x Xgenovo/assemble Xgenovo/finalize

ENV PATH="/NGStools/Xgenovo:$PATH"


