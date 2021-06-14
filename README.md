# polishCLR

## fetch pipeline

```
git clone https://github.com/isugifNF/polishCLR.git
cd polishCLR
```

## Install dependencies

For now, installing dependencies in a [miniconda](https://docs.conda.io/en/latest/miniconda.html) environment.

```
conda env create -f environment.yml -p ${PWD}/env/polishCLR_env
```

## Basic Run

**Ceres HPC**

```
module load nextflow
module load miniconda
source activate ${PWD}/env/polishCLR_env

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
executor >  slurm (1945)
[00/0890bc] process > bz_to_gz (1)            [100%] 1 of 1 ✔
[c0/97fd15] process > meryl_count_01 (1)      [100%] 2 of 2 ✔
[0c/aaa5c5] process > meryl_union_01          [100%] 1 of 1 ✔
[f2/12c259] process > MerquryQV_01 (1)        [100%] 1 of 1 ✔
[60/e66a84] process > pbmm2_index_02 (1)      [100%] 1 of 1 ✔
[6b/46feaa] process > pbmm2_align_02 (1)      [100%] 1 of 1 ✔
[10/03543b] process > create_windows_02 (1)   [100%] 1 of 1 ✔
[e4/d9a181] process > gcc_Arrow_02 (178)      [100%] 480 of 480 ✔
[20/c1adf5] process > merge_consensus_02      [100%] 1 of 1 ✔
[a1/6c5928] process > MerquryQV_03 (1)        [100%] 1 of 1 ✔
[6f/13e289] process > pbmm2_index_02b         [100%] 1 of 1 ✔
[75/18faa1] process > pbmm2_align_02b (1)     [100%] 1 of 1 ✔
[bf/324112] process > create_windows_02b      [100%] 1 of 1 ✔
[65/9cabdd] process > gcc_Arrow_02b (67)      [100%] 480 of 480 ✔
[a3/6925e8] process > merge_consensus_02b     [100%] 1 of 1 ✔
[99/ce7730] process > MerquryQV_03b (1)       [100%] 1 of 1 ✔
[d8/e99b27] process > align_shortreads_04 (1) [100%] 1 of 1 ✔
[6b/4f65ec] process > create_windows_04       [100%] 1 of 1 ✔
[10/5dfe3a] process > freebayes_04 (270)      [100%] 480 of 480 ✔
[c9/2478ab] process > combineVCF_04           [100%] 1 of 1 ✔
[94/a826ad] process > vcf_to_fasta_04 (1)     [100%] 1 of 1 ✔
[5c/4b6288] process > MerquryQV_05 (1)        [100%] 1 of 1 ✔
[f8/9fc043] process > align_shortreads_06 (1) [100%] 1 of 1 ✔
[e4/048537] process > create_windows_06 (1)   [100%] 1 of 1 ✔
[83/4a6b4f] process > freebayes_06 (444)      [100%] 480 of 480 ✔
[67/dda494] process > combineVCF_06           [100%] 1 of 1 ✔
[69/0f922a] process > vcf_to_fasta_06 (1)     [100%] 1 of 1 ✔
[f5/24cab8] process > MerquryQV_07 (1)        [100%] 1 of 1 ✔
Completed at: 12-Jun-2021 15:13:37
Duration    : 17h 32m 56s
CPU hours   : 82.1
Succeeded   : 1'945

```

[timeline.html](https://isugifnf.github.io/polishCLR/timeline.html) | [report.html](https://isugifnf.github.io/polishCLR/report.html)

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

PolishCLR_Results
|_ 00_Preprocess/               # Illumina bz2 reads converted to gz files
|_ 01_MerquryQV/                # quality value of primary assembly (before polishing)
|_ 02_ArrowPolish/              # polished with pacbio reads
|_ 03_MerquryQV/                # new quality value
|_ 02b_ArrowPolish/             # polished again with pacbio reads
|_ 03b_MerquryQV/               # new quality value
|_ 04_FreeBayesPolish/          # polished with illumina reads
|_ 05_MerquryQV/                 
|_ 06_FreeBayesPolish/    # <= should contain final assembly
|_ 07_MerquryQV/
|_ report.html
|_ timeline.html         # <= runtime for each step

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
