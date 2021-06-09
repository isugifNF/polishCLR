# polishCLR

## Install dependencies

For now, installing dependencies in a [miniconda](https://docs.conda.io/en/latest/miniconda.html) environment.

```
conda env create -f environment.yml -p ${PWD}/env/polishCLR_env
```

## Basic Run

**Ceres HPC**

```
module load nextflow
source activate ${PWD}/envs/polishCLR_env

nextflow run main.nf \
  --primary_assembly "/project/ag100pest/Pgos/RawData/3-unzip/all_p_ctg.fasta" \
  --illumina_reads "/project/ag100pest/Illumina_polishing/JAMU*{R1,R2}.fastq.bz2" \
  --pacbio_reads "/project/ag100pest/Pgos/RawData/m54334U_190823_194159.subreads.fasta" \
  --k "21" \
  -resume
```

pipeline for polishing CLRs

```
nextflow run main.nf -resume

N E X T F L O W  ~  version 20.10.0
Launching `main.nf` [mighty_franklin] - revision: 07a7c3d6ca
executor >  slurm (508)
[cd/518765] process > run_arrow (1) [100%] 508 of 508 ✔
Completed at: 19-May-2021 17:01:33
Duration    : 51m 51s
CPU hours   : 32.0
Succeeded   : 508
```
