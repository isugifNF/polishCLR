# polishCLR

## fetch pipeline

```
git clone https://github.com/isugifNF/polishCLR.git
cd polishCLR
```

<details><summary>print usage statement</summary>

  ```
  nextflow run main.nf --help
  
  N E X T F L O W  ~  version 21.04.0
Launching `main.nf` [condescending_gutenberg] - revision: 54771c06ba
  ----------------------------------------------------
                                \\---------//       
  ___  ___        _   ___  ___    \\-----//        
   |  (___  |  | / _   |   |_       \-//         
  _|_  ___) |__| \_/  _|_  |        // \        
                                  //-----\\       
                                //---------\\       
  isugifNF/polishCLR  v1.0.0       
----------------------------------------------------

Usage:
 The typical command for running the pipeline are as follows:
 nextflow run main.nf --primary_assembly "*fasta" --illumina_reads "*{1,2}.fastq.bz2" --pacbio_reads "*_subreads.bam" -resume

 Mandatory arguments:
 --primary_assembly             genome assembly fasta file to polish
 --illumina_reads               paired end illumina reads, to be used for Merqury QV scores, and freebayes polish primary assembly
 --pacbio_reads                 pacbio reads in bam format, to be used to arrow polish primary assembly
 --species                      if a string is given, rename the final assembly by species name [default:false]

 Optional modifiers
 --k                            kmer to use in MerquryQV scoring [default:21]
 --same_specimen                if illumina and pacbio reads are from the same specimin [default: true].
 --falcon_unzip                 if primary assembly has already undergone falcon unzip [default: false]. If true, will Arrow polish once instead of twice.

 Optional configuration arguments
 --parallel_app                 Link to parallel executable [default: 'parallel']
 --bzcat_app                    Link to bzcat executable [default: 'bzcat']
 --pigz_app                     Link to pigz executable [default: 'pigz']
 --meryl_app                    Link to meryl executable [default: 'meryl']
 --merqury_sh                   Link to merqury script [default: '$MERQURY/merqury.sh']
 --pbmm2_app                    Link to pbmm2 executable [default: 'pbmm2']
 --samtools_app                 Link to samtools executable [default: 'samtools']
 --gcpp_app                     Link to gcpp executable [default: 'gcpp']
 --bwamem2_app                  Link to bwamem2 executable [default: 'bwa-mem2']
 --freebayes_app                Link to freebayes executable [default: 'freebayes']
 --bcftools_app                 Link to bcftools executable [default: 'bcftools']
 --merfin_app                   Link to merfin executable [default: 'merfin']

 Optional arguments:
 --outdir                       Output directory to place final output [default: 'PolishCLR_Results']
 --clusterOptions               Cluster options for slurm or sge profiles [default slurm: '-N 1 -n 40 -t 04:00:00'; default sge: ' ']
 --threads                      Number of CPUs to use during each job [default: 40]
 --queueSize                    Maximum number of jobs to be queued [default: 50]
 --account                      Some HPCs require you supply an account name for tracking usage.  You can supply that here.
 --help                         This usage statement.
```
  
</details>

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

Merqury QV (quality value increases each polish, except 2nd arrow polish (investigate why)

```
=== original primary assembly
all_p_ctg    9064957 566335230	31.1459	0.000768078
=== Arrow polish
consensus 1601270	567139658	38.7086	0.000134629
=== Arrow polish
consensus 1638131	567090651	38.6093	0.000137745
=== Freebayes polish
consensus_01  1089557	567025929	40.3817	9.15852e-05
=== Freebayes polish
final_polished_assembly	1068195	567027516	40.4678	8.97877e-05
```

Final Polished Assembly QV = 40.4678
