FROM ubuntu:20.04

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -qqy \
    python3 \
    python3-pip \
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
    libssl-dev \
    unzip 


COPY ./liftover.py /usr/bin 
RUN chmod +x /usr/bin/liftover.py

RUN mkdir /resources
COPY ./hg19ToHg38.over.chain.gz /resources/


RUN pip3 install --upgrade pip
RUN pip3 install pandas==2.0.3
RUN pip3 install pyliftover==0.4.1

# Install htslib 1.20
WORKDIR /sources
RUN wget -O htslib-1.20.tar.bz2 https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2
RUN tar -xjvf htslib-1.20.tar.bz2
WORKDIR htslib-1.20/
RUN ./configure
RUN make
RUN make install


#RUN pip3 install trtools==6.0.2
# RUN git clone https://github.com/gymrek-lab/TRTools
# WORKDIR TRTools
# RUN pip3 install -e .

