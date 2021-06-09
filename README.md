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
  --pacbio_reads "/project/ag100pest/Pgos/RawData/m54334U_190823_194159.subreads.bam" \
  --k "21" \
  -resume
```

Current progress:

```
N E X T F L O W  ~  version 20.07.1
Launching `main.nf` [special_bartik] - revision: 7cb2c8ac31
executor >  slurm (10)
[53/50eca4] process > bz_to_gz (1)          [100%] 1 of 1, cached: 1 ✔
[23/3b4395] process > meryl_count_01 (1)    [100%] 2 of 2, cached: 2 ✔
[55/3e6914] process > meryl_union_01        [100%] 1 of 1, cached: 1 ✔
[a1/5ce12d] process > MerquryQV_01 (1)      [100%] 1 of 1, cached: 1 ✔
[fb/0af4a2] process > pbmm2_index_01 (1)    [100%] 1 of 1, cached: 1 ✔
[a6/1c9e97] process > pbmm2_align_01 (1)    [100%] 1 of 1, cached: 1 ✔
[53/4071bc] process > create_windows_01 (1) [100%] 1 of 1, cached: 1 ✔
[21/114bb7] process > gcc_Arrow_01 (12)     [  0%] 0 of 480
[-        ] process > merge_consensus_01    -
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
