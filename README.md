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
executor >  slurm (2)
[53/50eca4] process > bz_to_gz (1)            [100%] 1 of 1, cached: 1 ✔
[d7/24a868] process > meryl_count_01 (1)      [100%] 2 of 2, cached: 2 ✔
[55/3e6914] process > meryl_union_01          [100%] 1 of 1, cached: 1 ✔
[a1/5ce12d] process > MerquryQV_01 (1)        [100%] 1 of 1, cached: 1 ✔
[fb/0af4a2] process > pbmm2_index_01 (1)      [100%] 1 of 1, cached: 1 ✔
[a6/1c9e97] process > pbmm2_align_01 (1)      [100%] 1 of 1, cached: 1 ✔
[53/4071bc] process > create_windows_01 (1)   [100%] 1 of 1, cached: 1 ✔
[59/294ac0] process > gcc_Arrow_01 (478)      [100%] 480 of 480, cached: 480 ✔
[d8/b6c8d2] process > merge_consensus_01      [100%] 1 of 1, cached: 1 ✔
[13/1938bf] process > MerquryQV_02 (1)        [100%] 1 of 1, cached: 1 ✔
[9d/97a779] process > align_shortreads_01 (1) [100%] 1 of 1, cached: 1 ✔
[1d/e6da33] process > create_windows_02       [100%] 1 of 1, cached: 1 ✔
[99/1d2fcf] process > freebayes_01 (480)      [100%] 480 of 480, cached: 480 ✔
[d2/f8cd26] process > combineVCF_01           [100%] 1 of 1, cached: 1 ✔
[51/3f7ca7] process > vcf_to_fasta_01 (1)     [100%] 1 of 1, cached: 1 ✔
[56/3a97a1] process > MerquryQV_03 (1)        [100%] 1 of 1, cached: 1 ✔
[2d/86c139] process > align_shortreads_02 (1) [  0%] 0 of 1
[5d/94903c] process > create_windows_03 (1)   [100%] 1 of 1 ✔
[-        ] process > freebayes_02            -
[-        ] process > combineVCF_02           -
[-        ] process > vcf_to_fasta_02         -
[-        ] process > MerquryQV_04            -

```

Output Directory

```
PolishCLR_Results/
  |_ 00_Preprocess/              # illumina bz2 reads converted to gz files
  |_ 01_MerquryQV/
  |   |_ *.qv                    # quality value of original
  |_ 02_ArrowPolish/
  |   |_ consensus.fasta         # polished with pacbio reads
  |_ 03_MerquryQV/               # quality of arrow polished assembly
  |_ 04_FreeBayesPolish/
  |   |_ consensus_new.fasta     # Polished with illumina reads 
  |_ 05_MerquryQV/               # quality of freebayes polished assembly
  |_ 06_FreeBayesPolish/         # 2nd round ilumina polish
  |   |_ consensus_new.fasta     # <======= Final assemblly!
  |_ 07_MerquryQV/               # quality of final assembly

```

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