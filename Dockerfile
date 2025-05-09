FROM mambaorg/micromamba:2.1.1

ARG USER_ID
ARG GROUP_ID
ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN addgroup --gid $GROUP_ID user && \
    adduser --disabled-password --gecos '' --uid $USER_ID --gid $GROUP_ID user

COPY environment.yml /environment.yml
COPY environment_non_conda.sh /environment_non_conda.sh

RUN --mount=type=cache,target=/root/.cache/micromamba micromamba install -y -n base -f ./environment.yml -p /srv/onpdb_env && \
    micromamba clean --all --yes

ENV PATH=/srv/onpdb_env/bin:$PATH
ENV CONDA_DEFAULT_ENV=/srv/onpdb_env

WORKDIR /srv/onpdb

SHELL ["/bin/bash", "-c"]

RUN --mount=type=bind,source=environment_non_conda.sh,target=/tmp/environment_non_conda.sh micromamba run -p /srv/onpdb_env /tmp/environment_non_conda.sh

USER user

RUN micromamba init bash && \
    echo "micromamba activate /srv/onpdb_env" >> ~/.bashrc