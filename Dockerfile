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

RUN conda run -p /srv/onpdb_env /environment_non_conda.sh

USER user

RUN conda init bash
RUN echo "conda activate /srv/onpdb_env" >> ~/.bashrc
