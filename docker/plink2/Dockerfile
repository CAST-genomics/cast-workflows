FROM ubuntu:20.04

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -qqy \
    python3 \
    python3-pip \
    wget \
    build-essential \
    unzip \
    && rm -rf /var/lib/apt/lists/*


# Install Python dependencies
RUN pip install --no-cache-dir pandas

# Copy python script 
COPY convert_phenotype_plink.py /usr/bin/
RUN chmod +x /usr/bin/convert_phenotype_plink.py


# Download and install PLINK2
RUN wget https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20241020.zip \
    && unzip plink2_linux_x86_64_20241020.zip\
    && mv plink2 /usr/bin/ \
    && chmod +x /usr/bin/plink2 \
    && rm plink2_linux_x86_64_20241020.zip
    