FROM ubuntu:18.04

# Installation
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt-get install -y --no-install-recommends apt-utils
RUN apt-get upgrade -y
RUN apt-get install -y software-properties-common
RUN add-apt-repository universe
RUN apt-get update
RUN apt-get install -qqy wget
# Trying to install conda. Decided not to install.
#RUN wget https://repo.anaconda.com/archive/Anaconda2-2018.12-Linux-x86_64.sh
#RUN bash -y Anaconda2-2018.12-Linux-x86_64.sh
RUN apt-get install -qqy gcc
RUN apt-get install -qqy libbz2-dev
RUN apt-get install -qqy unzip
RUN apt-get install -qqy wget
RUN apt-get install -qqy curl
RUN apt-get install -qqy time
RUN apt-get install -qqy awscli
RUN apt-get install -qqy bzip2
RUN apt-get install -qqy python
RUN apt-get install -qqy libz-dev s3fs

# Install pip
RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py -o get-pip.py
RUN python get-pip.py


# Install samtools 1.13
# Install dependencies
RUN apt-get install -qqy make
RUN apt-get install -qqy libncurses-dev
RUN apt-get install -qqy liblzma-dev
RUN apt-get install -qqy autoconf libcurl4-openssl-dev
# Download, configure and install
RUN wget https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2
RUN tar -vxjf samtools-1.13.tar.bz2
WORKDIR samtools-1.13
RUN bash ./configure --enable-gcs --enable-libcurl
RUN make; make install
WORKDIR /

RUN apt-get install -qqy muscle python-tk

# Download advntr
# Download the latest updates on AdVNTR
RUN apt-get install -qqy git
RUN git clone -b aou_advntr https://github.com/mehrdadbakhtiari/adVNTR
WORKDIR adVNTR
RUN pwd

# Upgrade pip
#RUN python --version
#RUN pip --version
#RUN pip install --upgrade pip
#RUN pip --version
RUN apt-get install -qqy python-dev
RUN apt-get install -qqy build-essential

# Install pip dependencies
RUN pip install numpy==1.15.0
RUN pip install joblib==0.14.1
RUN pip install scipy==1.1.0
RUN pip install scikit-learn==0.16.1
RUN pip install Keras==2.3.1
RUN pip install cython==0.29.24
RUN pip install biopython==1.76
RUN pip install pysam==0.20.0
RUN pip install networkx==1.11
RUN pip install protobuf==3.17.3
RUN pip install tensorflow==1.13.2


# Add VNTR database and hg38 reference genome for seemless usage.
# hg38_VNTRs_by_TRF.db is a much larger database (826MB) compared to the p_vntrs_g_vntrs.db which is 19MB.
#COPY hg38_VNTRs_by_TRF.db vntr_db/hg38_VNTRs_by_TRF.db
COPY p_vntrs_g_vntrs.db vntr_db/p_vntrs_g_vntrs.db
COPY disease_associated_all_vntrs_v5_comma.txt vntr_db/phenotype_associated_vntrs_comma.txt
COPY extract_target_region.py /adVNTR

RUN ls -lh .
RUN ls -lh /adVNTR
RUN ls -lh /adVNTR/vntr_db/

# Install advntr
RUN make; make install
RUN python setup.py install

# Try running it to check version
RUN advntr genotype --help

