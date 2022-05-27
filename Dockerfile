FROM mambaorg/micromamba:0.19.1
USER root
COPY --chown=micromamba:micromamba environment.yml /tmp/env.yaml
RUN apt-get update
RUN apt-get install -y git curl python3-pip
RUN micromamba env create -y -f /tmp/env.yaml && \
    micromamba clean --all --yes
RUN ln -s /bin/micromamba /bin/conda

ARG MAMBA_DOCKERFILE_ACTIVATE=1 
RUN python -c 'import uuid; print(uuid.uuid4())' > /tmp/my_uuid

RUN echo "micromamba activate polishCLR_env" >> ~/.bashrc
RUN echo "export PATH=$PATH:${MAMBA_ROOT_PREFIX}/envs/polishCLR_env/bin" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]
ENV PATH "$PATH:${MAMBA_ROOT_PREFIX}/envs/polishCLR_env/bin"
