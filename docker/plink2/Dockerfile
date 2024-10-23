FROM ubuntu:20.04

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -qqy \
    python3 \
    python3-pip \
    wget \
    build-essential \
    curl \
    unzip \
    && rm -rf /var/lib/apt/lists/*


# Install Python dependencies
RUN pip3 install pandas

# Copy python script 
COPY convert_phenotype_plink.py /usr/bin/
RUN chmod +x /usr/bin/convert_phenotype_plink.py


# Download and install PLINK2
RUN wget https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20241020.zip \
    && unzip plink2_linux_x86_64_20241020.zip\
    && mv plink2 /usr/bin/ \
    && chmod +x /usr/bin/plink2 \
    && rm plink2_linux_x86_64_20241020.zip

# Download and install gcloud
RUN curl -sSL https://dl.google.com/dl/cloudsdk/release/google-cloud-sdk.tar.gz -o /tmp/google-cloud-sdk.tar.gz && \
    mkdir -p /usr/local/gcloud && \
    tar -C /usr/local/gcloud -xvf /tmp/google-cloud-sdk.tar.gz && \
    /usr/local/gcloud/google-cloud-sdk/install.sh --quiet && \
    rm -rf /tmp/google-cloud-sdk.tar.gz