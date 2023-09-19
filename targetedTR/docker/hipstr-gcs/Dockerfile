FROM ubuntu:20.04

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -qqy \
    awscli \
    zlib1g-dev libbz2-dev libz-dev liblzma-dev bzip2 \
    libcurl4-openssl-dev libssl-dev \
    python3 python3-pip \
    make gcc g++ \
    unzip \
    libfuse-dev fuse \
    git-all

RUN apt-get install -qqy \
    autoconf

# Install HipSTR
RUN git clone -b gcs https://github.com/gymrek-lab/HipSTR.git 
WORKDIR HipSTR/
RUN make version
WORKDIR lib/htslib
RUN autoreconf
RUN ./configure --enable-gcs --enable-s3 --enable-libcurl --without-libdeflate
WORKDIR ../../
RUN make
RUN cp HipSTR /usr/local/bin/
WORKDIR ..