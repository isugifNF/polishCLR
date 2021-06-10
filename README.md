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
executor >  slurm (482)
[53/50eca4] process > bz_to_gz (1)            [100%] 1 of 1, cached: 1 ✔
[d7/24a868] process > meryl_count_01 (2)      [100%] 2 of 2, cached: 2 ✔
[55/3e6914] process > meryl_union_01          [100%] 1 of 1, cached: 1 ✔
[a1/5ce12d] process > MerquryQV_01 (1)        [100%] 1 of 1, cached: 1 ✔
[fb/0af4a2] process > pbmm2_index_01 (1)      [100%] 1 of 1, cached: 1 ✔
[a6/1c9e97] process > pbmm2_align_01 (1)      [100%] 1 of 1, cached: 1 ✔
[53/4071bc] process > create_windows_01 (1)   [100%] 1 of 1, cached: 1 ✔
[8d/598162] process > gcc_Arrow_01 (477)      [100%] 480 of 480, cached: 480 ✔
[d8/b6c8d2] process > merge_consensus_01      [100%] 1 of 1, cached: 1 ✔
[13/1938bf] process > MerquryQV_02 (1)        [100%] 1 of 1, cached: 1 ✔
[9d/97a779] process > align_shortreads_01 (1) [100%] 1 of 1 ✔
[1d/e6da33] process > create_windows_02       [100%] 1 of 1, cached: 1 ✔
[16/f31f64] process > freebayes_01 (309)      [100%] 480 of 480 ✔
[cf/17e899] process > combineVCF_01           [100%] 1 of 1 ✔
Completed at: 10-Jun-2021 17:28:55
Duration    : 1h 16m 37s
CPU hours   : 46.2 (91.8% cached)
Succeeded   : 482
Cached      : 491
```

Continue from `consensus.vcf` either later today or tomorrow.

<details><summary>prior run</summary>

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

</details>