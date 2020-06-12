FROM continuumio/anaconda3
COPY environment_loose.yml /environment_loose.yml
RUN conda-env create -f=./environment_loose.yml  -p /srv/onpdb_env