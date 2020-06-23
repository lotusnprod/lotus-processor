FROM continuumio/miniconda3

ARG USER_ID
ARG GROUP_ID
RUN addgroup --gid $GROUP_ID user
RUN adduser --disabled-password --gecos '' --uid $USER_ID --gid $GROUP_ID user

COPY environment_loose.yml /environment_loose.yml
COPY environment_non_conda.sh /environment_non_conda.sh

RUN conda-env create -f=./environment_loose.yml  -p /srv/onpdb_env
ENV PATH /srv/onpdb_env/bin:$PATH
ENV CONDA_DEFAULT_ENV /srv/onpdb_env

RUN mkdir /srv/onpdb
WORKDIR /srv/onpdb

SHELL ["/bin/bash", "-c"]

RUN cd /tmp && curl -L https://github.com/gnames/gnfinder/releases/download/v0.11.1/gnfinder-v0.11.1-linux.tar.gz | tar xz && mv gnfinder /usr/local/bin
RUN conda run -p /srv/onpdb_env /environment_non_conda.sh

USER user
RUN echo "conda activate /srv/onpdb_env" >> ~/.bashrc