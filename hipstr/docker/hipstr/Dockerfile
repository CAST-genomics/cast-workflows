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

# Install HipSTR
RUN git clone -b v1.0.0 https://github.com/gymrek-lab/HipSTR.git 
WORKDIR HipSTR/
RUN make version
RUN make
RUN cp HipSTR /usr/local/bin/
WORKDIR ..

CMD ["/bin/bash"]

