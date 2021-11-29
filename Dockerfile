FROM mambaorg/micromamba

# === Update installers
RUN apt-get update && apt-get install -y git && apt-get -y install mafft
RUN apt-get install -y curl && apt-get install -y python3-pip

# === Pull repo
RUN git clone https://github.com/isugifNF/polishCLR.git

# === install any conda packages
RUN cd polishCLR; mamba env create -f environment.yml

# === merfin not on conda
RUN cd; git clone https://github.com/arangrhie/merfin.git
RUN cd merfin/src
RUN make -j 8
RUN ./merfin --version

# === install nextflow?
RUN curl -s https://get.nextflow.io | bash
RUN mv nextflow ~/bin/.

# Rough outline, still need to build to check if this works