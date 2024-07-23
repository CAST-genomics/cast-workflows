FROM ubuntu:20.04

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -qqy \
    python3 \
    python3-pip \
    git-all wget \
    curl \
    make \
    cmake \
    gcc \
    autoconf \
    automake \
    libtool \
    subversion \
    pkg-config \
    openjdk-21-jdk 


#default-jre \


# Download and install gcloud
RUN curl https://dl.google.com/dl/cloudsdk/release/google-cloud-sdk.tar.gz > /tmp/google-cloud-sdk.tar.gz
RUN mkdir -p /usr/local/gcloud \
  && tar -C /usr/local/gcloud -xvf /tmp/google-cloud-sdk.tar.gz \
  && /usr/local/gcloud/google-cloud-sdk/install.sh
ENV PATH $PATH:/usr/local/gcloud/google-cloud-sdk/bin


RUN wget "https://faculty.washington.edu/browning/beagle/beagle.27May24.118.jar" -O beagle.jar




