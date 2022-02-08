## Nextflow polishCLR pipeline
*polishCLR* is a [nextflow](https://www.nextflow.io/) workflow for polishing genome assemblies (improving accuary) generated with noisy PacBio reads using accurate, short Illumina reads. It implements the best practices descbribed by the Vertebrate Genome Project (VGP) Assembly community (Rhie et al. 2021) and extends these for use-cases we found common in the [Ag100Pest Genome Initiative](http://i5k.github.io/ag100pest).

The polishCLR workflow can be easily initiated from three input cases:
- Case 1: An unresolved primary assembly with associated contigs (the output of FALCON 2-asm) or without (e.g., the output of Canu or wtdbg2). 
- Case 2: A haplotype-resolved but unpolished set (e.g., the output of FALCON-Unzip 3-unzip). 
- **Case 3: IDEAL! A haplotype-resolved, CLR long-read, Arrow-polished set of primary and alternate contigs (e.g., the output of FALCON-Unzip 4-polish).** 

We strongly reccomend including the organellular genome to improve the polishing of nuclear mitochondrial or plasmid pseudogenes (Howe et al., 2021). Organelle genomes should be generated and polished separately for best results. You could use the mitochondrial companion to polishCLR, [polishCLRmt](https://github.com/Ag100Pest/Ag100MitoPolishCLR) or [mitoVGP](https://github.com/gf777/mitoVGP) (Formenti et al., 2021). 

To allow for the inclusion of scaffolding before final polishing  and increase the potential for gap-filling across correctly oriented scaffolded contigs, the core workflow is divided into two steps, controlled by a `--step` parameter flag. 

<p align="center">
	<img src="https://github.com/isugifNF/polishCLR/blob/main/Figure01.svg" width="750" />
  
</p>

You can view a more complete vizualization of the pipleine in [Supp. Fig. S1](https://github.com/isugifNF/polishCLR/blob/main/FigureS01.svg)

## Documentation
You can find more details on the usage below. These also include a simple [step-by-step] tutorial to run the analyses on your own genomes.

## Table of Contents

- [Installation](#Installation)
- [Basic Run](#Basic-Run)
- [Outputs](#Outputs)
- [Trio Input](#Trio-Input)
- [References](#References)


## Installation
Fetch pipeline

```
git clone https://github.com/isugifNF/polishCLR.git
cd polishCLR
```
Install dependencies in a [miniconda](https://docs.conda.io/en/latest/miniconda.html) environment.

```
conda env create -f environment.yml -p ${PWD}/env/polishCLR_env
```

## Basic Run

Print usage statement

  ```
N E X T F L O W  ~  version 21.04.2
Launching `main.nf` [lonely_liskov] - revision: 6a81970115
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
   --illumina_reads               paired end illumina reads, to be used for Merqury QV scores, and freebayes polish primary assembly
   --pacbio_reads                 pacbio reads in bam format, to be used to arrow polish primary assembly
   --mitochondrial_assembly       mitocondrial assembly will be concatinated to the assemblies before polishing [default: false]

   Either FALCON (or FALCON Unzip) assembly:
   --primary_assembly             genome assembly fasta file to polish
   --alternate_assembly           if alternate/haplotig assembly file is provided, will be concatinated to the primary assembly before polishing [default: false]
   --falcon_unzip                 if primary assembly has already undergone falcon unzip [default: false]. If true, will Arrow polish once instead of twice.

   Or TrioCanu assembly
   --paternal_assembly            paternal genome assembly fasta file to polish
   --maternal_assembly            maternal genome assembly fasta file to polish

   Pick Step 1 (arrow, purgedups) or Step 2 (arrow, freebayes, freebayes)
   --step                         Run step 1 or step 2 (default: 1)

   Optional modifiers
   --species                      if a string is given, rename the final assembly by species name [default:false]
   --k                            kmer to use in MerquryQV scoring [default:21]
   --same_specimen                if illumina and pacbio reads are from the same specimin [default: true].
   --meryldb                      path to a prebuilt meryl database, built from the illumina reads. If not provided, tehen build.

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


**Example**

```
source activate ${PWD}/env/polishCLR_env

SP_DIR=/project/ag100pest/Helicoverpa_zea_male

## === Step 1 (up to purge_dups)
nextflow run /project/ag100pest/software/polishCLR/main.nf \
  --primary_assembly "${SP_DIR}/Falcon/4-polish/cns-output/cns_p_ctg.fasta" \
  --alternate_assembly "${SP_DIR}/Falcon/4-polish/cns-output/cns_h_ctg.fasta" \
  --mitochondrial_assembly "${SP_DIR}/MT_Contig/Hzea_mtDNA_contig.fasta" \
  --illumina_reads "${SP_DIR}/PolishingData/JDRP*{R1,R2}.fastq.bz2" \
  --pacbio_reads "${SP_DIR}/RawData/m*.subreads.bam" \
  --species "Hzea" \
  --k "21" \
  --falcon_unzip true \
  --step 1 \
  --busco_lineage "lepideoptera_odb10" \
  -resume \
  -profile ceres

### I often submit these as two separate slurm scripts

## === Step 2 
 nextflow run /project/ag100pest/software/polishCLR/main.nf \
  --primary_assembly "PolishCLR_Results/Step_1/02_Purge_Dups/primary_purged.fa" \
  --alternate_assembly "PolishCLR_Results/Step_1/02_Purge_Dups/haps_purged.fa" \
  --mitochondrial_assembly "${SP_DIR}/MT_Contig/Hzea_mtDNA_contig.fasta" \
  --illumina_reads "${SP_DIR}/PolishingData/JDRP*{R1,R2}.fastq.bz2" \
  --pacbio_reads "${SP_DIR}/RawData/m*.subreads.bam" \
  --species "Hzea" \
  --k "21" \
  --falcon_unzip true \
  --step 2 \
  -resume \
  -profile ceres
```
### Outputs
Key Ouputs are found in
```
- PolishCLR_Results/Step_1/02_BUSCO/primary_purged/ ## BUSCO stats after purge_dups
- PolishCLR_Results/Step_1/02_Purge_Dups/primary_purged.fa ## Primary assembly after Step 1
- PolishCLR_Results/Step_1/02_Purge_Dups/haps_purged.fa ## Alternate assembly after Step 1
- PolishCLR_Results/Step_2/*_bbstat/ ## Length and number distributions through each process of Step 2
- PolishCLR_Results/Step_2/*_QV/ ## QV and Completeness stats with merfin and merqury for each process of Step 2
- PolishCLR_Results/Step_2/06_FreeBayesPolish/p_Step_2_06_FreeBayesPolish_consensus.fasta ## Final primary fasta (no mito)
- PolishCLR_Results/Step_2/06_FreeBayesPolish/a_Step_2_06_FreeBayesPolish_consensus.fasta ## Final alternate fasta (no mito)
```

### Example Output Logging
**Step 1**
```
[ea/f78cb8] process > RENAME_PRIMARY (1) [100%] 1 of 1, cached: 1 ✔
[53/770510] process > MERGE_FILE_00 (1)  [100%] 1 of 1, cached: 1 ✔
[34/a4b7ac] process > bz_to_gz (1)       [100%] 1 of 1, cached: 1 ✔
[3e/0af766] process > meryl_count (1)    [100%] 2 of 2, cached: 2 ✔
[76/a26404] process > meryl_union        [100%] 1 of 1, cached: 1 ✔
[7c/d4916c] process > meryl_peak         [100%] 1 of 1, cached: 1 ✔
[ef/c58401] process > MerquryQV_00 (1)   [100%] 1 of 1, cached: 1 ✔
[cd/3ad9c7] process > bbstat_00 (1)      [100%] 1 of 1, cached: 1 ✔
[3e/dd32f5] process > bam_to_fasta (1)   [100%] 1 of 1, cached: 1 ✔
[00/1b1177] process > SPLIT_FILE_02 (1)  [100%] 1 of 1, cached: 1 ✔
[78/df2d30] process > PURGE_DUPS_02 (1)  [100%] 1 of 1, cached: 1 ✔
[5f/dec926] process > BUSCO (2)          [100%] 2 of 2 ✔

Completed at: 21-Dec-2021 10:45:52
Duration    : 51m 32s
CPU hours   : 3.9 (63.3% cached)
Succeeded   : 2
Cached      : 12
```

**Step 2** 
```
executor >  local (1), slurm (2681)
[d9/b48b2c] process > RENAME_PRIMARY (1)             [100%] 1 of 1 ✔
[25/b9acec] process > MERGE_FILE_00 (1)              [100%] 1 of 1 ✔
[94/ff9db2] process > bz_to_gz (1)                   [100%] 1 of 1 ✔
[2c/bd2fc0] process > meryl_count (1)                [100%] 2 of 2 ✔
[00/b4d5aa] process > meryl_union                    [100%] 1 of 1 ✔
[62/116067] process > meryl_peak                     [100%] 1 of 1 ✔
[76/7ecb59] process > MerquryQV_00 (1)               [100%] 1 of 1 ✔
[12/ee4395] process > bbstat_00 (1)                  [100%] 1 of 1 ✔
[93/7d56c3] process > ARROW_04:create_windows (1)    [100%] 1 of 1 ✔
[3a/15bd43] process > ARROW_04:pbmm2_index (1)       [100%] 1 of 1 ✔
[05/75aae9] process > ARROW_04:pbmm2_align (1)       [100%] 1 of 1 ✔
[6b/92ea05] process > ARROW_04:gcpp_arrow (3)        [100%] 882 of 882 ✔
[55/29df60] process > ARROW_04:meryl_genome (1)      [100%] 1 of 1 ✔
[0f/8eb0eb] process > ARROW_04:combineVCF (1)        [100%] 1 of 1 ✔
[9b/a8433d] process > ARROW_04:reshape_arrow (1)     [100%] 1 of 1 ✔
[2f/bb106e] process > ARROW_04:merfin_polish (1)     [100%] 1 of 1 ✔
[e5/609325] process > ARROW_04:vcf_to_fasta (1)      [100%] 1 of 1 ✔
[d6/0b30fa] process > MerquryQV_04 (1)               [100%] 1 of 1 ✔
[4c/dd5ce5] process > bbstat_04 (1)                  [100%] 1 of 1 ✔
[5c/7d00a5] process > FREEBAYES_05:create_windows... [100%] 1 of 1 ✔
[62/269cb4] process > FREEBAYES_05:meryl_genome (1)  [100%] 1 of 1 ✔
[43/329d5d] process > FREEBAYES_05:align_shortrea... [100%] 1 of 1 ✔
[15/ece2b2] process > FREEBAYES_05:freebayes (754)   [100%] 882 of 882 ✔
[ab/47faff] process > FREEBAYES_05:combineVCF (1)    [100%] 1 of 1 ✔
[8b/5461c2] process > FREEBAYES_05:merfin_polish (1) [100%] 1 of 1 ✔
[ea/0a6fb3] process > FREEBAYES_05:vcf_to_fasta (1)  [100%] 1 of 1 ✔
[4a/c67e0a] process > MerquryQV_05 (1)               [100%] 1 of 1 ✔
[57/f661d4] process > bbstat_05 (1)                  [100%] 1 of 1 ✔
[a5/5e8166] process > FREEBAYES_06:create_windows... [100%] 1 of 1 ✔
[ea/ccf140] process > FREEBAYES_06:meryl_genome (1)  [100%] 1 of 1 ✔
[40/d1cbd4] process > FREEBAYES_06:align_shortrea... [100%] 1 of 1 ✔
[f3/a7920f] process > FREEBAYES_06:freebayes (514)   [100%] 882 of 882 ✔
[ff/c3c29b] process > FREEBAYES_06:combineVCF (1)    [100%] 1 of 1 ✔
[9f/462b6c] process > FREEBAYES_06:merfin_polish (1) [100%] 1 of 1 ✔
[a0/3218e0] process > FREEBAYES_06:vcf_to_fasta (1)  [100%] 1 of 1 ✔
[f8/0f5b7c] process > MerquryQV_06 (1)               [100%] 1 of 1 ✔
[e8/09d9db] process > bbstat_06 (1)                  [100%] 1 of 1 ✔
[12/04a0db] process > SPLIT_FILE_07 (1)              [100%] 1 of 1 ✔
Completed at: 23-Dec-2021 11:45:16
Duration    : 21h 33m 40s
CPU hours   : 195.0
Succeeded   : 2'682

```


[timeline.html](https://isugifnf.github.io/polishCLR/timeline_ceres.html) | [report.html](https://isugifnf.github.io/polishCLR/report_ceres.html)

We also developed a pipeline to work with Trio data, but this is less tested
<details><summary>Trio Usage</summary>
### Trio Input
**TrioCanu Assembly**

```
nextflow run main.nf \
  -stub-run \
  -profile local \
  --paternal_assembly "data/pri.fasta" \
  --maternal_assembly "data/alt.fasta" \
  --mitochondrial_assembly "data/mit.fasta" \
  --illumina_reads "data/readname_{R1,R2}.fastq.bz" \
  --pacbio_reads "data/*subreads.fasta" \
  --outdir "TrioPolish_Polish"
```

TrioCanu Step 1 Output

```
N E X T F L O W  ~  version 21.04.3
Launching `main.nf` [elegant_wiles] - revision: bfeb7f0055
executor >  local (34)
[02/be6f35] process > MERGE_FILE_TRIO (1)           [100%] 2 of 2 ✔
[a0/5c5042] process > meryl_count (1)               [100%] 2 of 2 ✔
[f0/9c536a] process > meryl_union                   [100%] 1 of 1 ✔
[9b/ce96e3] process > meryl_peak                    [100%] 1 of 1 ✔
[e4/1c3238] process > MerquryQV_01 (2)              [100%] 2 of 2 ✔
[36/b1494c] process > bbstat_01 (2)                 [100%] 2 of 2 ✔
[4c/46aa68] process > ARROW_02:create_windows (1)   [100%] 1 of 1 ✔
[c0/e9bb62] process > ARROW_02:pbmm2_index (1)      [100%] 1 of 1 ✔
[44/0e1bce] process > ARROW_02:pbmm2_align (1)      [100%] 1 of 1 ✔
[d3/e97fcd] process > ARROW_02:gcc_Arrow (3)        [100%] 3 of 3 ✔
[90/d306b3] process > ARROW_02:merge_consensus (1)  [100%] 1 of 1 ✔
[d5/13127b] process > ARROW_02b:create_windows (1)  [100%] 1 of 1 ✔
[92/302898] process > ARROW_02b:pbmm2_index (1)     [100%] 1 of 1 ✔
[a7/1a0dc2] process > ARROW_02b:pbmm2_align (1)     [100%] 1 of 1 ✔
[4e/8994ef] process > ARROW_02b:gcc_Arrow (3)       [100%] 3 of 3 ✔
[a3/3fc599] process > ARROW_02b:merge_consensus (1) [100%] 1 of 1 ✔
[1a/83e513] process > MerquryQV_03 (2)              [100%] 2 of 2 ✔
[cc/399faa] process > bbstat_03 (2)                 [100%] 2 of 2 ✔
[ca/734f45] process > SPLIT_FILE_03p                [100%] 1 of 1 ✔
[40/61755d] process > PURGE_DUPS_TRIOp (1)          [100%] 1 of 1 ✔
[1a/cd0464] process > BUSCO (1)                     [100%] 1 of 1 ✔
[6e/797242] process > SPLIT_FILE_03m                [100%] 1 of 1 ✔
[39/f44778] process > PURGE_DUPS_TRIOm (1)          [100%] 1 of 1 ✔
[c5/2fa18b] process > BUSCO_mat (1)                 [100%] 1 of 1 ✔
```

**TrioCanu Assembly**

```
nextflow run main.nf \
  -stub-run \
  -profile local \
  --paternal_assembly "data/pri.fasta" \
  --maternal_assembly "data/alt.fasta" \
  --mitochondrial_assembly "data/mit.fasta" \
  --illumina_reads "data/readname_{R1,R2}.fastq.bz" \
  --pacbio_reads "data/*subreads.fasta" \
  --outdir "TrioPolish_Polish" \
  --steptwo true
```

```
N E X T F L O W  ~  version 21.04.3
Launching `main.nf` [gloomy_lichterman] - revision: bfeb7f0055
executor >  local (82)
[d7/056083] process > MERGE_FILE_TRIO (1)                [100%] 2 of 2 ✔
[93/51edce] process > meryl_count (1)                    [100%] 2 of 2 ✔
[ed/5527db] process > meryl_union                        [100%] 1 of 1 ✔
[e6/e804bc] process > meryl_peak                         [100%] 1 of 1 ✔
[63/c47356] process > MerquryQV_01 (2)                   [100%] 2 of 2 ✔
[c1/bb67de] process > bbstat_01 (2)                      [100%] 2 of 2 ✔
[8b/490e0f] process > ARROW_04:create_windows (1)        [100%] 1 of 1 ✔
[a3/a5e0b1] process > ARROW_04:pbmm2_index (1)           [100%] 1 of 1 ✔
[ed/a2259b] process > ARROW_04:pbmm2_align (1)           [100%] 1 of 1 ✔
[a3/020fa4] process > ARROW_04:gcc_Arrow (2)             [100%] 3 of 3 ✔
[eb/631b37] process > ARROW_04:meryl_genome (1)          [100%] 1 of 1 ✔
[01/ebc79b] process > ARROW_04:combineVCF (1)            [100%] 1 of 1 ✔
[c5/84d378] process > ARROW_04:reshape_arrow (1)         [100%] 1 of 1 ✔
[69/b0cde0] process > ARROW_04:merfin_polish (1)         [100%] 1 of 1 ✔
[f3/8b6850] process > ARROW_04:vcf_to_fasta (1)          [100%] 1 of 1 ✔
[14/effa19] process > ARROW_04b:create_windows (1)       [100%] 1 of 1 ✔
[49/f3c8fb] process > ARROW_04b:pbmm2_index (1)          [100%] 1 of 1 ✔
[dc/093b9d] process > ARROW_04b:pbmm2_align (1)          [100%] 1 of 1 ✔
[15/a58989] process > ARROW_04b:gcc_Arrow (3)            [100%] 3 of 3 ✔
[d5/d4b0dd] process > ARROW_04b:meryl_genome (1)         [100%] 1 of 1 ✔
[cb/77a45b] process > ARROW_04b:combineVCF (1)           [100%] 1 of 1 ✔
[e1/984069] process > ARROW_04b:reshape_arrow (1)        [100%] 1 of 1 ✔
[de/7d54db] process > ARROW_04b:merfin_polish (1)        [100%] 1 of 1 ✔
[20/de921d] process > ARROW_04b:vcf_to_fasta (1)         [100%] 1 of 1 ✔
[af/8fb965] process > MerquryQV_05 (2)                   [100%] 2 of 2 ✔
[57/17e453] process > bbstat_05 (2)                      [100%] 2 of 2 ✔
[b3/db7637] process > FREEBAYES_06:create_windows (1)    [100%] 1 of 1 ✔
[8f/561e8e] process > FREEBAYES_06:meryl_genome (1)      [100%] 1 of 1 ✔
[fa/737fd4] process > FREEBAYES_06:align_shortreads (1)  [100%] 1 of 1 ✔
[4c/a9817c] process > FREEBAYES_06:freebayes (2)         [100%] 3 of 3 ✔
[79/9b75f0] process > FREEBAYES_06:combineVCF (1)        [100%] 1 of 1 ✔
[6a/aa4244] process > FREEBAYES_06:merfin_polish (1)     [100%] 1 of 1 ✔
[b6/fb2986] process > FREEBAYES_06:vcf_to_fasta (1)      [100%] 1 of 1 ✔
[26/272e2a] process > FREEBAYES_06b:create_windows (1)   [100%] 1 of 1 ✔
[ae/28f8fb] process > FREEBAYES_06b:meryl_genome (1)     [100%] 1 of 1 ✔
[ae/068a62] process > FREEBAYES_06b:align_shortreads (1) [100%] 1 of 1 ✔
[e2/8884d6] process > FREEBAYES_06b:freebayes (3)        [100%] 3 of 3 ✔
[bd/97676a] process > FREEBAYES_06b:combineVCF (1)       [100%] 1 of 1 ✔
[8e/6b69ad] process > FREEBAYES_06b:merfin_polish (1)    [100%] 1 of 1 ✔
[1f/9dbbb2] process > FREEBAYES_06b:vcf_to_fasta (1)     [100%] 1 of 1 ✔
[5e/626d17] process > MerquryQV_07 (2)                   [100%] 2 of 2 ✔
[61/e076cd] process > bbstat_07 (2)                      [100%] 2 of 2 ✔
[7c/f18e5f] process > FREEBAYES_08:create_windows (1)    [100%] 1 of 1 ✔
[3a/effebc] process > FREEBAYES_08:meryl_genome (1)      [100%] 1 of 1 ✔
[f3/ae1183] process > FREEBAYES_08:align_shortreads (1)  [100%] 1 of 1 ✔
[f5/365d05] process > FREEBAYES_08:freebayes (1)         [100%] 3 of 3 ✔
[b7/064817] process > FREEBAYES_08:combineVCF (1)        [100%] 1 of 1 ✔
[ac/c326bb] process > FREEBAYES_08:merfin_polish (1)     [100%] 1 of 1 ✔
[cd/8fba43] process > FREEBAYES_08:vcf_to_fasta (1)      [100%] 1 of 1 ✔
[ee/928768] process > FREEBAYES_08b:create_windows (1)   [100%] 1 of 1 ✔
[f4/386ef4] process > FREEBAYES_08b:meryl_genome (1)     [100%] 1 of 1 ✔
[b5/e7edeb] process > FREEBAYES_08b:align_shortreads (1) [100%] 1 of 1 ✔
[50/ae5389] process > FREEBAYES_08b:freebayes (1)        [100%] 3 of 3 ✔
[27/e057bb] process > FREEBAYES_08b:combineVCF (1)       [100%] 1 of 1 ✔
[59/6387f0] process > FREEBAYES_08b:merfin_polish (1)    [100%] 1 of 1 ✔
[56/539a81] process > FREEBAYES_08b:vcf_to_fasta (1)     [100%] 1 of 1 ✔
[97/dffa16] process > MerquryQV_09 (2)                   [100%] 2 of 2 ✔
[1b/a81a7a] process > bbstat_09 (2)                      [100%] 2 of 2 ✔
[83/592eb3] process > SPLIT_FILE_09p                     [100%] 1 of 1 ✔
[b4/aa4ddb] process > SPLIT_FILE_09m                     [100%] 1 of 1 ✔
```
</details>

# References
Rhie, A., McCarthy, S. A., Fedrigo, O., Damas, J., Formenti, G., Koren, S., Uliano-Silva, M., Chow, W., Fungtammasan, A., Kim, J., Lee, C., Ko, B. J., Chaisson, M., Gedman, G. L., Cantin, L. J., Thibaud-Nissen, F., Haggerty, L., Bista, I., Smith, M., . . . Jarvis, E. D. (2021). Towards complete and error-free genome assemblies of all vertebrate species. Nature, 592(7856), 737-746. https://doi.org/10.1038/s41586-021-03451-0 

Rhie, A., Walenz, B. P., Koren, S., & Phillippy, A. M. (2020). Merqury: reference-free quality, completeness, and phasing assessment for genome assemblies. Genome Biol, 21(1), 245. https://doi.org/10.1186/s13059-020-02134-9.

Formenti, G., Rhie, A., Balacco, J., Haase, B., Mountcastle, J., Fedrigo, O., Brown, S., Capodiferro, M. R., Al-Ajli, F. O., Ambrosini, R., Houde, P., Koren, S., Oliver, K., Smith, M., Skelton, J., Betteridge, E., Dolucan, J., Corton, C., Bista, I., . . . Vertebrate Genomes Project, C. (2021). Complete vertebrate mitogenomes reveal widespread repeats and gene duplications. Genome Biol, 22(1), 120. https://doi.org/10.1186/s13059-021-02336-9 
