FROM mambaorg/micromamba:0.23.3
USER root
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/env.yaml
RUN sed -i 's/polishCLR_env/base/g' /tmp/env.yaml
RUN apt-get update
RUN apt-get install -y git curl python3-pip
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes
ENV PATH "$PATH:$MAMBA_ROOT_PREFIX/envs/polishCLR_env/bin"