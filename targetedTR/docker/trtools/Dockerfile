FROM ubuntu:20.04

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -qqy \
    python3 \
    python3-pip

RUN pip3 install --upgrade pip
RUN pip3 install pandas
RUN git clone  -b conf_ref  https://github.com/gymreklab/TRTools
WORKDIR TRTools/
RUN pip install --upgrade pip
RUN pip install -e .






