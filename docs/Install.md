# Installation

Fetch pipeline

```
git clone https://github.com/isugifNF/polishCLR.git
cd polishCLR
```

For ag100pest projects, we have been installing dependencies in a [miniconda](https://docs.conda.io/en/latest/miniconda.html) environment. Dependencies may also be installed via [docker](https://www.docker.com/) or [singularity](https://sylabs.io/guides/3.0/user-guide/quick_start.html).

### Miniconda

Install dependencies in a [miniconda](https://docs.conda.io/en/latest/miniconda.html) environment.

```
[[ -d env ]] || mkdir env
conda env create -f environment.yml -p ${PWD}/env/polishCLR_env
```

### Docker

Start up [docker](https://docs.docker.com/get-docker/) and pull the [csiva2022/polishclr:latest](https://hub.docker.com/r/csiva2022/polishclr) image.

```
docker pull csiva2022/polishclr:latest
```

Run the polishCLR pipeline with an added `-with-docker csiva2022/polishclr:latest ` parameter. See [nextflow docker run documentation](https://www.nextflow.io/docs/latest/docker.html#how-it-works) for more information.

### Singularity

Install dependencies as a [singularity](https://sylabs.io/guides/3.0/user-guide/quick_start.html) image.

```
singularity pull polishclr.sif docker://csiva2022/polishclr:latest
```

Run the polishCLR pipeline with an added `-with-singularity polishclr.sif ` parameter. See [nextflow singularity run documentation](https://www.nextflow.io/docs/latest/singularity.html#how-it-works) for more information.
