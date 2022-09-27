FROM mambaorg/micromamba:0.25.0
USER root
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/env.yaml
RUN sed -i 's/polishCLR_env/base/g' /tmp/env.yaml
RUN apt-get update
RUN apt-get install -y git curl python3-pip procps
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes
RUN rm /opt/conda/lib/jvm/bin/java
RUN ln -s /opt/conda/pkgs/openjdk-11.0.15-h1e1ecb3_2/bin/java /opt/conda/lib/jvm/bin/java
ENV PATH "$PATH:$MAMBA_ROOT_PREFIX/bin"
