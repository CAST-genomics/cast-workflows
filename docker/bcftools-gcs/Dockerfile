FROM ubuntu:20.04

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -qqy \
    python3 \
    python3-pip \
    curl \
    build-essential \
    gcc \
    wget \
    autoconf \
    pkg-config \
    libbz2-dev \
    liblzma-dev \
    libncurses5 libncurses5-dev \
    libcurl4-gnutls-dev \
    libnss3-dev \
    libssl-dev 


# Make sources directory
RUN mkdir /sources

# Download and install gcloud
RUN curl https://dl.google.com/dl/cloudsdk/release/google-cloud-sdk.tar.gz > /tmp/google-cloud-sdk.tar.gz
RUN mkdir -p /usr/local/gcloud \
  && tar -C /usr/local/gcloud -xvf /tmp/google-cloud-sdk.tar.gz \
  && /usr/local/gcloud/google-cloud-sdk/install.sh
ENV PATH $PATH:/usr/local/gcloud/google-cloud-sdk/bin


# Install htslib 1.20
WORKDIR /sources
RUN wget -O htslib-1.20.tar.bz2 https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2
RUN tar -xjvf htslib-1.20.tar.bz2
WORKDIR htslib-1.20/
RUN ./configure
RUN make
RUN make install


#Install bcftools 1.20
WORKDIR /sources
RUN wget -O bcftools-1.20.tar.bz2 https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2
RUN tar -xjvf bcftools-1.20.tar.bz2
WORKDIR bcftools-1.20
RUN ./configure --enable-gcs --enable-libcurl --enable-s3 
RUN make install


