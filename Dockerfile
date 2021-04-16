FROM continuumio/miniconda3

ARG USER_ID
ARG GROUP_ID
RUN addgroup --gid $GROUP_ID user
RUN adduser --disabled-password --gecos '' --uid $USER_ID --gid $GROUP_ID user

COPY environment.yml /environment.yml
COPY environment_non_conda.sh /environment_non_conda.sh

RUN conda-env create -f=./environment.yml  -p /srv/onpdb_env
ENV PATH /srv/onpdb_env/bin:$PATH
ENV CONDA_DEFAULT_ENV /srv/onpdb_env

RUN mkdir /srv/onpdb
WORKDIR /srv/onpdb

SHELL ["/bin/bash", "-c"]

## AR: I made commands available in the makefile to replace following do not know how to include them here

RUN cd /tmp && curl -L https://github.com/gnames/gnfinder/releases/download/v0.11.1/gnfinder-v0.11.1-linux.tar.gz | tar xz && mv gnfinder /usr/local/bin
RUN cd /tmp && curl -L https://github.com/gnames/gnverify/releases/download/v0.1.0/gnverify-v0.1.0-linux.tar.gz | tar xz && mv gnverify /usr/local/bin

## TODO


RUN conda run -p /srv/onpdb_env /environment_non_conda.sh

USER user

RUN conda init bash
RUN echo "conda activate /srv/onpdb_env" >> ~/.bashrc