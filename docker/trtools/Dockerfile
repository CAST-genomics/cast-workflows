FROM ubuntu:20.04

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -qqy \
    python3 \
    python3-pip \
    git-all wget 


RUN pip3 install --upgrade pip
RUN pip3 install pandas
RUN pip3 install trtools==6.1.0
#RUN git clone https://github.com/gymrek-lab/TRTools
#WORKDIR TRTools
#RUN pip3 install -e .

