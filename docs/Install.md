# Installation

This pipeline will require [nextflow](https://www.nextflow.io/docs/latest/getstarted.html). The rest fo the dependencies can be installed via miniconda, Docker, or Singularity

```
nextflow -version

nextflow run isugifNF/polishCLR -r main --help
```

For ag100pest projects, we have been installing dependencies in a [miniconda](https://docs.conda.io/en/latest/miniconda.html) environment. Dependencies may also be installed via [docker](https://www.docker.com/) or [singularity](https://sylabs.io/guides/3.0/user-guide/quick_start.html).

Since nextflow caches pulled pipelines in a `~/.nextflow` folder somewhere on your system, to update the pipeline we recommend:

```
nextflow drop isugifNF/polishCLR
nextflow pull isugifNF/polishCLR -r main
```

### Miniconda

Install dependencies in a [miniconda](https://docs.conda.io/en/latest/miniconda.html) environment.

```
wget https://raw.githubusercontent.com/isugifNF/polishCLR/main/environment.yml

[[ -d env ]] || mkdir env
conda env create -f environment.yml -p ${PWD}/env/polishCLR_env

conda activate env/polishCLR_env

nextflow run isugifNF/polishCLR -r main \
  --check_software
```

### Docker

Start up [docker](https://docs.docker.com/get-docker/) and pull the [csiva2022/polishclr:latest](https://hub.docker.com/r/csiva2022/polishclr) image.

```
docker pull csiva2022/polishclr:latest

# Option 1
docker run -it csiva2022/polishclr:latest \
  nextflow run isugifNF/polishCLR -r main \
  --check_software

# Option 2
nextflow run isugifNF/polishCLR -r main \
  --check_software \
  -profile docker

# Option 3
nextflow run isugifNF/polishCLR -r main \
  --check_software \
  --with-docker csiva2022/polishclr:latest
```

Run the polishCLR pipeline with an added `-with-docker csiva2022/polishclr:latest ` parameter. See [nextflow docker run documentation](https://www.nextflow.io/docs/latest/docker.html#how-it-works) for more information.

### Singularity

Install dependencies as a [singularity](https://sylabs.io/guides/3.0/user-guide/quick_start.html) image.

```
singularity pull polishclr.sif docker://csiva2022/polishclr:latest

# Option 1:
singularity exec polishclr.sif \
  nextflow run isugifNF/polishCLR -r main \
  --check_software

# Option 2:
nextflow run isugifNF/polishCLR -r main \
  --check_software \
  -profile singularity

# Option 3:
nextflow run isugifNF/polishCLR -r main \
  --check_software \
  -with-singularity polishclr.sif
```

Run the polishCLR pipeline with an added `-with-singularity polishclr.sif ` parameter. See [nextflow singularity run documentation](https://www.nextflow.io/docs/latest/singularity.html#how-it-works) for more information.

