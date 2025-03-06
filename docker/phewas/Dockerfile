FROM ubuntu:20.04

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -qqy \
    python3 \
    python3-pip \
    curl \
    git

# Clone the directory
RUN git clone https://github.com/CAST-genomics/cast-workflows.git
WORKDIR cast-workflows/phewas
RUN git checkout sara_main_2
RUN git pull origin sara_main_2
RUN pip install cyvcf2 pandas

# Install PheTK
RUN pip install PheTK --upgrade

# To resolve missing required cpu features error.
RUN pip install polars-lts-cpu

# Download and install gcloud
RUN curl https://dl.google.com/dl/cloudsdk/release/google-cloud-sdk.tar.gz > /tmp/google-cloud-sdk.tar.gz
RUN mkdir -p /usr/local/gcloud \
  && tar -C /usr/local/gcloud -xvf /tmp/google-cloud-sdk.tar.gz \
    && /usr/local/gcloud/google-cloud-sdk/install.sh
ENV PATH $PATH:/usr/local/gcloud/google-cloud-sdk/bin
