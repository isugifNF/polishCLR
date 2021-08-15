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

<details><summary>Install `merfin` if not already in a module</summary>
  
https://github.com/arangrhie/merfin

```
alloc -N 1 -n 8 -p scavenger -t 04:00:00
module load git
module load gcc/8.1.0
git clone https://github.com/arangrhie/merfin.git
cd merfin/src
make -j 8
./merfin --version
```
  
</details>

<details><summary>Install `genomescope2.0` if not already in a module</summary>
 
https://github.com/tbenavi1/genomescope2.0

```
alloc -N 1 -n 8 -p scavenger -t 04:00:00
git clone https://github.com/tbenavi1/genomescope2.0.git
cd genomescope2.0
emacs install.R       # edit to set 'local_lib_path = "/project/ag100pest/software/R_libs/"'
module load r/3.4.1
Rscript install.R
emacs genomescope.R   # add the lib path to each library 'library(XXXXX, lib.loc = local_lib_path)'
Rscript genomescope.R --help
```
  
</details>

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
executor >  slurm (967)
[45/e5018c] process > bz_to_gz (1)                 [100%] 1 of 1, cached: 1 ✔
[a7/a5b5a1] process > meryl_count (2)              [100%] 2 of 2, cached: 2 ✔
[c2/81d8d5] process > meryl_union                  [100%] 1 of 1, cached: 1 ✔
[d3/984a04] process > MerquryQV_01 (1)             [100%] 1 of 1, cached: 1 ✔
[40/aa96c3] process > bbstat_01 (1)                [100%] 1 of 1, cached: 1 ✔
[03/9b95b5] process > ARROW_02:create_windows (1)  [100%] 1 of 1, cached: 1 ✔
[62/7a60a8] process > ARROW_02:pbmm2_index (1)     [100%] 1 of 1, cached: 1 ✔
[5c/5f0217] process > ARROW_02:pbmm2_align (1)     [100%] 1 of 1 ✔
[1c/3b58a6] process > ARROW_02:gcc_Arrow (95)      [100%] 480 of 480 ✔
[67/97354d] process > ARROW_02:merge_consensus (1) [100%] 1 of 1 ✔
[e7/b6a605] process > MerquryQV_03 (1)             [100%] 1 of 1 ✔
[41/b87aea] process > bbstat_03 (1)                [100%] 1 of 1 ✔
[b2/2dc1e1] process > ARROW_04:create_windows (1)  [100%] 1 of 1 ✔
[eb/f11dc2] process > ARROW_04:pbmm2_index (1)     [100%] 1 of 1 ✔
[eb/9c76e8] process > ARROW_04:pbmm2_align (1)     [100%] 1 of 1 ✔
[ea/6f9cb5] process > ARROW_04:gcc_Arrow (134)     [100%] 480 of 480 ✔ // pause, add in merfin as check
[-        ] process > ARROW_04:merge_consensus      -
[-        ] process > MerquryQV_05                  -
[-        ] process > bbstat_05                     -
[-        ] process > FREEBAYES_06:create_windows   -
[-        ] process > FREEBAYES_06:align_shortreads -
[-        ] process > FREEBAYES_06:freebayes        -
[-        ] process > FREEBAYES_06:combineVCF       -
[-        ] process > FREEBAYES_06:vcf_to_fasta     -
[-        ] process > MerquryQV_07                  -
[-        ] process > bbstat_07                     -
[-        ] process > FREEBAYES_08:create_windows   -
[-        ] process > FREEBAYES_08:align_shortreads -
[-        ] process > FREEBAYES_08:freebayes        -
[-        ] process > FREEBAYES_08:combineVCF       -
[-        ] process > FREEBAYES_08:vcf_to_fasta     -
[-        ] process > MerquryQV_09                  -
[-        ] process > bbstat_09                     -
```

[timeline.html](https://isugifnf.github.io/polishCLR/timeline.html) | [report.html](https://isugifnf.github.io/polishCLR/report.html)

Output Directory

```
PolishCLR_Results/
  |_ 00_Preprocess/               # Illumina bz2 reads converted to gz files
  |_ 01_QV/                       # quality value of primary assembly (before polishing)
  |  |_ MerquryQV/                # Merqury quality value and histogram plots
  |  |_ bbstat/                   # bbstat quality value
  |_ 02_ArrowPolish/              # polished with pacbio reads
  |  |_ gccruns/                  # subfolder of vcf and fasta files per window
  |  |_ 2_consensus.fasta         # arrow polished new consensus sequence
  |_ 03_QV/                       # New merqury and bbstat quality values in subfolders
  |_ 04_ArrowPolish/              # polished again with pacbio reads
  |_ 05_QV/
  |_ 06_FreeBayesPolish/          # polished with illumina reads
  |_ 07_QV/
  |_ 08_FreeBayesPolish/          # <= should contain final polished assembly 8_consensus.fasta
  |_ 09_QV/ 
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
