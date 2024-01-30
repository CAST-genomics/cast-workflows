FROM ubuntu:18.04

# Installation
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install -qqy gcc
RUN apt-get install -qqy libbz2-dev
RUN apt-get install -qqy unzip
RUN apt-get install -qqy wget
RUN apt-get install -qqy awscli
RUN apt-get install -qqy bzip2
RUN apt-get install -qqy python python-pip
RUN apt-get install -qqy libz-dev s3fs
RUN apt-get install -qqy samtools muscle python-tk
RUN apt-get install -qqy time

# Download advntr
#RUN wget https://github.com/mehrdadbakhtiari/adVNTR/archive/1.5.0.tar.gz
#RUN tar -xzvf 1.5.0.tar.gz
# Download the latest updates on AdVNTR
RUN apt-get install -qqy git
RUN git clone https://github.com/mehrdadbakhtiari/adVNTR
WORKDIR adVNTR
RUN pwd

# Install pip dependencies
RUN pip install numpy==1.15.0
RUN pip install joblib==0.14.1
RUN pip install scipy==1.1.0
RUN pip install scikit-learn==0.16.1
RUN pip install Keras==2.3.1
RUN pip install cython==0.29.24
RUN pip install biopython==1.76
RUN pip install pysam==0.9.1.4
RUN pip install networkx==1.11
RUN pip install protobuf==3.17.3
RUN pip install tensorflow==1.13.2


# Add VNTR database and hg38 reference genome for seemless usage.
COPY hg38_VNTRs_by_TRF.db vntr_db/hg38_VNTRs_by_TRF.db
COPY disease_associated_all_vntrs_v4_comma.txt vntr_db/phenotype_associated_vntrs_comma.txt

RUN ls -lh .
RUN ls -lh /adVNTR
RUN ls -lh /adVNTR/vntr_db/

# Install advntr
RUN make; make install
RUN python setup.py install

# Try running it to check version
RUN advntr genotype --help

