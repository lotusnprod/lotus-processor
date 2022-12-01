FROM mambaorg/micromamba:latest

ARG USER_ID
ARG GROUP_ID
ARG MAMBA_DOCKERFILE_ACTIVATE=1  # otherwise python will not be found
RUN addgroup --gid $GROUP_ID user
RUN adduser --disabled-password --gecos '' --uid $USER_ID --gid $GROUP_ID user

COPY environment.yml /environment.yml
COPY environment_non_conda.sh /environment_non_conda.sh

RUN micromamba install -y -n base -f ./environment.yml -p /srv/onpdb_env && \
    micromamba clean --all --yes 
ENV PATH /srv/onpdb_env/bin:$PATH
ENV CONDA_DEFAULT_ENV /srv/onpdb_env

RUN mkdir /srv/onpdb
WORKDIR /srv/onpdb

SHELL ["/bin/bash", "-c"]

RUN micromamba run -p /srv/onpdb_env /environment_non_conda.sh

USER user

RUN micromamba init bash
RUN echo "micromamba activate /srv/onpdb_env" >> ~/.bashrc
