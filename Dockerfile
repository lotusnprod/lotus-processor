FROM continuumio/miniconda3
COPY environment_loose.yml /environment_loose.yml
RUN conda-env create -f=./environment_loose.yml  -p /srv/onpdb_env
ENV PATH /srv/onpdb_env/bin:$PATH
ENV CONDA_DEFAULT_ENV onpdb_env
RUN mkdir /srv/onpdb
WORKDIR /srv/onpdb
RUN cd /tmp && curl -L https://github.com/gnames/gnfinder/releases/download/v0.11.1/gnfinder-v0.11.1-linux.tar.gz | tar xz && mv gnfinder /usr/local/bin
