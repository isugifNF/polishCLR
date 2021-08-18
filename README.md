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
module load merfin

SP_DIR=/project/ag100pest/Pectinophora_gossypiella_male/

nextflow run main.nf \
  --primary_assembly "$SP_DIR/RawData/3-unzip/all_p_ctg.fasta" \
  --mito_assembly "$SP_DIR/MT_Contig/Pgos/Pgos_MitoFinder_mitfi_Final_Results/Pgos_mtDNA_contig.fasta" \
  --illumina_reads "$SP_DIR/Illumina_polishing/JAMU*{R1,R2}.fastq.bz2" \
  --pacbio_reads "$SP_DIR/RawData/m54334U_190823_194159.subreads.bam" \
  --k "21" \
  -resume
```

Current progress:

```
N E X T F L O W  ~  version 20.07.1
Launching `main.nf` [special_bartik] - revision: 7cb2c8ac31
executor >  slurm (1954)
[4b/d5fb33] process > addMito (1)                    [100%] 1 of 1, cached: 1 ✔
[9d/d3d19f] process > bz_to_gz (1)                   [100%] 1 of 1, cached: 1 ✔
[de/ef0633] process > meryl_count (2)                [100%] 2 of 2, cached: 2 ✔
[fb/7ef659] process > meryl_union                    [100%] 1 of 1, cached: 1 ✔
[23/ae01da] process > meryl_peak                     [100%] 1 of 1, cached: 1 ✔
[bb/96e7da] process > MerquryQV_01 (1)               [100%] 1 of 1, cached: 1 ✔  // QV=31.1461
[ae/c81f3c] process > bbstat_01 (1)                  [100%] 1 of 1, cached: 1 ✔
[1c/f9f4a2] process > ARROW_02:create_windows (1)    [100%] 1 of 1, cached: 1 ✔
[f6/a11eb2] process > ARROW_02:pbmm2_index (1)       [100%] 1 of 1, cached: 1 ✔
[99/55ebfc] process > ARROW_02:pbmm2_align (1)       [100%] 1 of 1 ✔
[8c/e62875] process > ARROW_02:gcc_Arrow (107)       [100%] 481 of 481 ✔
[0e/e86ee7] process > ARROW_02:merge_consensus (1)   [100%] 1 of 1 ✔
[a7/2a2420] process > MerquryQV_03 (1)               [100%] 1 of 1 ✔   // QV=38.7089
[3f/8c66e2] process > bbstat_03 (1)                  [100%] 1 of 1 ✔
[05/220810] process > ARROW_04:create_windows (1)    [100%] 1 of 1 ✔
[ab/1cab8c] process > ARROW_04:pbmm2_index (1)       [100%] 1 of 1 ✔
[0c/1a92d7] process > ARROW_04:pbmm2_align (1)       [100%] 1 of 1 ✔
[4e/80a70b] process > ARROW_04:gcc_Arrow (68)        [100%] 481 of 481 ✔
[1a/530d23] process > ARROW_04:meryl_genome (1)      [100%] 1 of 1 ✔
[79/d30118] process > ARROW_04:combineVCF_arrow (1)  [100%] 1 of 1 ✔
[5a/ca838c] process > ARROW_04:reshape_arrow (1)     [100%] 1 of 1 ✔
[b7/2cc664] process > ARROW_04:merfin_polish_arro... [100%] 1 of 1 ✔
[25/c8d4cd] process > ARROW_04:vcf_to_fasta_arrow... [100%] 1 of 1 ✔  // QV=41.0985
[71/21f8ea] process > MerquryQV_05 (1)               [100%] 1 of 1 ✔
[17/8bee1a] process > bbstat_05 (1)                  [100%] 1 of 1 ✔
[99/dd5c4f] process > FREEBAYES_06:create_windows... [100%] 1 of 1 ✔
[ea/62b6f0] process > FREEBAYES_06:meryl_genome_f... [100%] 1 of 1 ✔
[0d/b74cde] process > FREEBAYES_06:align_shortrea... [100%] 1 of 1 ✔
[39/2f8bbc] process > FREEBAYES_06:freebayes (335)   [100%] 481 of 481 ✔
[23/75ce2a] process > FREEBAYES_06:combineVCF (1)    [100%] 1 of 1 ✔
[e0/298ae8] process > FREEBAYES_06:merfin_polish (1) [100%] 1 of 1 ✔
[41/27460a] process > FREEBAYES_06:vcf_to_fasta (1)  [100%] 1 of 1 ✔
[48/cf8399] process > MerquryQV_07 (1)               [100%] 1 of 1 ✔  // QV=45.0154
[3e/5a07de] process > bbstat_07 (1)                  [100%] 1 of 1 ✔
[2b/515801] process > FREEBAYES_08:create_windows... [100%] 1 of 1 ✔
[23/53bd22] process > FREEBAYES_08:meryl_genome_f... [100%] 1 of 1 ✔
[db/2833ab] process > FREEBAYES_08:align_shortrea... [100%] 1 of 1 ✔
[35/71633c] process > FREEBAYES_08:freebayes (68)    [100%] 481 of 481 ✔
[57/6d06bb] process > FREEBAYES_08:combineVCF (1)    [100%] 1 of 1 ✔
[40/cc4a82] process > FREEBAYES_08:merfin_polish (1) [100%] 1 of 1 ✔
[f5/de1f18] process > FREEBAYES_08:vcf_to_fasta (1)  [100%] 1 of 1 ✔
[3c/fdab12] process > MerquryQV_09 (1)               [100%] 1 of 1 ✔  // QV=45.0468
[d3/009a91] process > bbstat_09 (1)                  [100%] 1 of 1 ✔
Completed at: 17-Aug-2021 23:14:48
Duration    : 17h 33m 38s
CPU hours   : 89.2 (0.9% cached)
Succeeded   : 1'954
Cached      : 10
```

<!--
```
executor >  slurm (967)
[45/e5018c] process > bz_to_gz (1)                 [100%] 1 of 1, cached: 1 ✔
[a7/a5b5a1] process > meryl_count (2)              [100%] 2 of 2, cached: 2 ✔
[c2/81d8d5] process > meryl_union                  [100%] 1 of 1, cached: 1 ✔
[d3/984a04] process > MerquryQV_01 (1)             [100%] 1 of 1, cached: 1 ✔ // QV score of 31.1459
[40/aa96c3] process > bbstat_01 (1)                [100%] 1 of 1, cached: 1 ✔
[03/9b95b5] process > ARROW_02:create_windows (1)  [100%] 1 of 1, cached: 1 ✔
[62/7a60a8] process > ARROW_02:pbmm2_index (1)     [100%] 1 of 1, cached: 1 ✔
[5c/5f0217] process > ARROW_02:pbmm2_align (1)     [100%] 1 of 1 ✔
[1c/3b58a6] process > ARROW_02:gcc_Arrow (95)      [100%] 480 of 480 ✔
[67/97354d] process > ARROW_02:merge_consensus (1) [100%] 1 of 1 ✔
[e7/b6a605] process > MerquryQV_03 (1)             [100%] 1 of 1 ✔  // QV score of 38.7086
[41/b87aea] process > bbstat_03 (1)                [100%] 1 of 1 ✔
[b2/2dc1e1] process > ARROW_04:create_windows (1)  [100%] 1 of 1 ✔
[eb/f11dc2] process > ARROW_04:pbmm2_index (1)     [100%] 1 of 1 ✔
[eb/9c76e8] process > ARROW_04:pbmm2_align (1)     [100%] 1 of 1 ✔
[ea/6f9cb5] process > ARROW_04:gcc_Arrow (134)     [100%] 480 of 480 ✔ 
[cf/53361f] process > ARROW_04:meryl_peak (1)      [100%] 1 of 1, cached: 1 ✔
[20/5f153d] process > ARROW_04:meryl_genome (1)    [100%] 1 of 1, cached: 1 ✔
[3e/7bef6a] process > ARROW_04:combineVCF_arrow (1)[100%] 1 of 1, cached: 1 ✔
[36/40b176] process > ARROW_04:reshape_arrow (1)   [100%] 1 of 1, cached: 1 ✔
[55/aca211] process > ARROW_04:merfin_polish_arro..[100%] 1 of 1 ✔
[d3/919202] process > ARROW_04:vcf_to_fasta_arrow..[100%] 1 of 1 ✔
[de/5a06a8] process > MerquryQV_05 (1)             [100%] 1 of 1 ✔   // Wow, QV score of 41.1038, use 2nd arrow polish and pass to FreeBayes
[a7/abfafd] process > bbstat_05 (1)                [100%] 1 of 1 ✔
[8a/17ed90] process > FREEBAYES_06:create_windows..[100%] 1 of 1, cached: 1 ✔
[4d/e8cd59] process > FREEBAYES_06:meryl_genome_f..[100%] 1 of 1, cached: 1 ✔
[20/1cbb5c] process > FREEBAYES_06:align_shortrea..[100%] 1 of 1, cached: 1 ✔
[74/6a11c5] process > FREEBAYES_06:freebayes (475) [100%] 480 of 480, cache...
[72/fb44a3] process > FREEBAYES_06:combineVCF (1)  [100%] 1 of 1, cached: 1 ✔
[5c/05ecb1] process > FREEBAYES_06:merfin_polish (1[100%] 1 of 1, cached: 1 ✔
[38/f8c31b] process > FREEBAYES_06:vcf_to_fasta (1)[100%] 1 of 1 ✔
[a3/c96d68] process > MerquryQV_07 (1)             [100%] 1 of 1 ✔  // QV without merfin 42.7133, QV with merfin 45.0122
[69/48a675] process > bbstat_07 (1)                [100%] 1 of 1 ✔
[60/91de55] process > FREEBAYES_08:create_windows..[100%] 1 of 1 ✔
[ee/927c53] process > FREEBAYES_08:meryl_genome_f..[100%] 1 of 1 ✔
[b1/a1dfef] process > FREEBAYES_08:align_shortrea... [100%] 1 of 1 ✔
[2f/ae39b7] process > FREEBAYES_08:freebayes (303)   [100%] 480 of 480 ✔
[12/a14839] process > FREEBAYES_08:combineVCF (1)    [100%] 1 of 1 ✔
[c0/8a3aeb] process > FREEBAYES_08:merfin_polish (1) [100%] 1 of 1 ✔
[40/0f9772] process > FREEBAYES_08:vcf_to_fasta (1)  [100%] 1 of 1 ✔
[9b/e157d4] process > MerquryQV_09 (1)               [100%] 1 of 1 ✔ // QV 45.0425
[d6/476547] process > bbstat_09 (1)                  [100%] 1 of 1 ✔
Completed at: 16-Aug-2021 23:26:44
Duration    : 3h 43m 34s
CPU hours   : 96.5 (93.9% cached)
Succeeded   : 491
Cached      : 1'468
```
-->

[timeline.html](https://isugifnf.github.io/polishCLR/timeline.html) | [report.html](https://isugifnf.github.io/polishCLR/report.html)

Output Directory

```
PolishCLR_Results/
  |_ 00_Preprocess/               # Illumina bz2 reads converted to gz files
  |_ 01_QV/                       # quality value of primary assembly (before polishing)
  |  |_ MerquryQV/                # Merqury quality value and histogram plots
  |  |_ bbstat/                   # bbstat quality value
  |  |_ merqury.qv                # <= text file with only the merqury qv value, can use as a checkpoint
  |_ 02_ArrowPolish/              # polished with pacbio reads
  |  |_ gccruns/                  # subfolder of vcf and fasta files per window
  |  |_ 2_consensus.fasta         # arrow polished new consensus sequence
  |_ 03_QV/                       # New merqury and bbstat quality values in subfolders. Continue with the assembly with the higher merqury.qv value
  |_ 04_ArrowPolish/              # polished again with pacbio reads
  |  |_ merfin/                   # merfin filtering if illumina and pacbio from same sample
  |_ 05_QV/
  |_ 06_FreeBayesPolish/          # polished with illumina reads
  |  |_ merfin/                   # merfin filtering
  |_ 07_QV/
  |_ 08_FreeBayesPolish/          # <= should contain final polished assembly 8_consensus.fasta
  |_ 09_QV/
  |_ report.html
  |_ timeline.html         # <= runtime for each step
```

<!--
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
-->

Final Polished Assembly QV = 45.0425
