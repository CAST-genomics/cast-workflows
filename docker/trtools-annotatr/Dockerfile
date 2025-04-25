FROM ubuntu:22.04

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -qqy \
    python3 \
    python3-pip \
    git-all wget 


RUN pip3 install --upgrade pip
RUN pip3 install pandas
#RUN pip3 install trtools==6.0.2
# RUN git clone https://github.com/gymrek-lab/TRTools
# WORKDIR TRTools
# RUN pip3 install -e .

RUN pip install --upgrade --force-reinstall git+https://github.com/yli091230/TRTools.git@master

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -qqy \
    curl \
    build-essential \
    gcc \
    autoconf \
    pkg-config \
    libbz2-dev \
    liblzma-dev \
    libncurses5 libncurses5-dev \
    libcurl4-gnutls-dev \
    libnss3-dev \
    libssl-dev \
    unzip

# Install htslib 1.20
WORKDIR /sources
RUN wget -O htslib-1.20.tar.bz2 https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2
RUN tar -xjvf htslib-1.20.tar.bz2
WORKDIR htslib-1.20/
RUN ./configure
RUN make
RUN make install
