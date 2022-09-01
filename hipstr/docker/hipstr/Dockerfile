FROM ubuntu:20.04

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -qqy \
    awscli \
    zlib1g-dev libbz2-dev libz-dev liblzma-dev bzip2 \
    libcurl4-openssl-dev libssl-dev \
    python3 python3-pip \
    make gcc g++ \
    unzip \
    libfuse-dev fuse \
    curl \
    wget

# Install HipSTR
#RUN wget  https://github.com/gymrek-lab/HipSTR/archive/refs/heads/main.zip -O HipSTR.zip
#RUN unzip  HipSTR.zip 
#WORKDIR HipSTR-main/
#RUN make
#RUN cp HipSTR /usr/local/bin/
#WORKDIR ..

# Install HipSTR
RUN wget https://github.com/tfwillems/HipSTR/archive/v0.7.tar.gz
RUN tar -xzvf v0.7.tar.gz
WORKDIR HipSTR-0.7/
RUN make
RUN cp HipSTR /usr/local/bin/
WORKDIR ..




CMD ["/bin/bash"]

