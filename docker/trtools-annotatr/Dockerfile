FROM continuumio/miniconda3


RUN apt-get update && apt-get install -y git

RUN git clone https://github.com/yli091230/TRTools 
WORKDIR TRTools
RUN conda env create -n trtools -f dev-env.yml
ENV PATH /opt/conda/envs/trtools/bin:$PATH
RUN conda run -n trtools poetry install
# Set the default command to bash to keep the container running
CMD ["/bin/bash"]

