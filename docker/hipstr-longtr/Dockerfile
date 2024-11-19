FROM ubuntu:20.04

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -qqy \
    autoconf \
    awscli \
    curl \
    zlib1g-dev libbz2-dev libz-dev liblzma-dev bzip2 \
    libcurl4-openssl-dev libssl-dev \
    python3 python3-pip \
    make gcc g++ \
    unzip \
    libfuse-dev fuse \
    git-all wget \
    cmake

# Download and install gcloud
RUN curl https://dl.google.com/dl/cloudsdk/release/google-cloud-sdk.tar.gz > /tmp/google-cloud-sdk.tar.gz
RUN mkdir -p /usr/local/gcloud \
  && tar -C /usr/local/gcloud -xvf /tmp/google-cloud-sdk.tar.gz \
  && /usr/local/gcloud/google-cloud-sdk/install.sh
ENV PATH $PATH:/usr/local/gcloud/google-cloud-sdk/bin

# Download HTSLIB
RUN wget https://github.com/samtools/htslib/releases/download/1.19/htslib-1.19.tar.bz2
RUN tar -xvf htslib-1.19.tar.bz2
WORKDIR htslib-1.19
RUN ./configure --enable-gcs --enable-s3 --enable-libcurl
RUN make -j && make install

# Install HipSTR GCS htslib update version
RUN git clone https://github.com/gymrek-lab/HipSTR.git 
WORKDIR HipSTR/
RUN git checkout 2c83bed 
RUN make version
RUN make -j
RUN cp HipSTR /usr/local/bin/
WORKDIR ..
RUN ldconfig

# Install LongTR
RUN git clone https://github.com/gymrek-lab/longtr
WORKDIR longtr/lib
RUN git clone https://github.com/rvaser/spoa
WORKDIR ..

RUN make
RUN cp LongTR /usr/local/bin/
