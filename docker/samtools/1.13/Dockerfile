FROM ubuntu:18.04

# Installation
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt-get install -y --no-install-recommends apt-utils
RUN apt-get upgrade -y
RUN apt-get install -qqy wget \
                          gcc \
                          libbz2-dev \
                          bzip2 \
                          curl \
                          python \
                          libz-dev \
                          s3fs

# Install pip
RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py -o get-pip.py
RUN python get-pip.py

# Install samtools 1.13
# Install dependencies
RUN apt-get install -qqy make \
                          libncurses-dev \
                          liblzma-dev \
                          autoconf \
                          libcurl4-openssl-dev
# Download, configure and install
RUN wget https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2
RUN tar -vxjf samtools-1.13.tar.bz2
WORKDIR samtools-1.13
RUN bash ./configure --enable-gcs --enable-libcurl
RUN make; make install
WORKDIR /

RUN apt-get install -qqy muscle python-tk

