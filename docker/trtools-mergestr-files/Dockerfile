FROM ubuntu:20.04

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -qqy \
    git-all \
    python3 \
    python3-pip

RUN pip3 install --upgrade pip
RUN git clone https://github.com/gymrek-lab/trtools.git 
WORKDIR trtools
RUN git checkout f07a99d
RUN pip install -e .

# Install dependencies
RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get -qq install -y \
    autoconf \
    build-essential \
    libbz2-dev \
    liblzma-dev \
    libncurses5 libncurses5-dev \
    pkg-config \
    wget \
    zlib1g-dev

# Make sources directory
RUN mkdir /sources

# Install VCFTools
WORKDIR /sources
RUN wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz
RUN tar -xzvf vcftools-0.1.16.tar.gz
WORKDIR vcftools-0.1.16
RUN ./autogen.sh && ./configure && make && make install
