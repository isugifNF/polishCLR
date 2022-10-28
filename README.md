## Nextflow polishCLR pipeline

[![DOI](https://zenodo.org/badge/375112950.svg)](https://zenodo.org/badge/latestdoi/375112950) <a href="https://hub.docker.com/r/csiva2022/polishclr">
      <img alt="isugifnf PolishCLR version" src="https://img.shields.io/docker/v/csiva2022/polishclr?label=%F0%9F%90%8B%20%20%20docker%3Apolishclr">
  </a> [![Build Status](https://github.com/isugifNF/polishCLR/actions/workflows/stubtest.yml/badge.svg?branch=main)](https://github.com/isugifNF/polishCLR/actions/workflows/stubtest.yml)

*polishCLR* is a [nextflow](https://www.nextflow.io/) workflow for polishing genome assemblies (improving accuracy) generated with noisy PacBio reads using accurate, short Illumina reads. It implements the best practices described by the Vertebrate Genome Project (VGP) Assembly community (Rhie et al. 2021) and extends these for use-cases we found common in the [Ag100Pest Genome Initiative](http://i5k.github.io/ag100pest). This workflow was developed as part of the USDA-ARS Ag100Pest Initiative. The authors thank members of the USDA-ARS Ag100Pest Team and SCINet Virtual Resource Support Core (VRSC) for fruitful discussions and troubleshooting throughout the development of this workflow. 

The polishCLR workflow can be easily initiated from three input cases:
- Case 1: An unresolved primary assembly with associated contigs (the output of FALCON 2-asm) or without (e.g., the output of Canu or wtdbg2). 
- Case 2: A haplotype-resolved but unpolished set (e.g., the output of FALCON-Unzip 3-unzip). 
- **Case 3: IDEAL! A haplotype-resolved, CLR long-read, Arrow-polished set of primary and alternate contigs (e.g., the output of FALCON-Unzip 4-polish).** 

We strongly recommend including the organellular genome to improve the polishing of nuclear mitochondrial or plasmid pseudogenes (Howe et al., 2021). Organelle genomes should be generated and polished separately for best results. You could use the mitochondrial companion to polishCLR, [polishCLRmt](https://github.com/Ag100Pest/Ag100MitoPolishCLR) or [mitoVGP](https://github.com/gf777/mitoVGP) (Formenti et al., 2021). 

To allow for the inclusion of scaffolding before final polishing  and increase the potential for gap-filling across correctly oriented scaffolded contigs, the core workflow is divided into two steps, controlled by a `--step` parameter flag. 

<p align="center">
	<img src="https://github.com/isugifNF/polishCLR/blob/main/docs/imgs/Figure01.svg" width="750" />
</p>

You can view a more complete visualization of the pipeline in [Supp. Fig. S1](https://github.com/isugifNF/polishCLR/blob/main/FigureS01.svg)

## Documentation
You can find more details on the usage below. These also include a simple [basic run](#Basic-Run) to run the analyses on your own genomes.

## Table of Contents

- [Installation](#Installation)
- [Basic Run](#Basic-Run)
- [Outputs](#Outputs)
- [Trio Input](#Trio-Input)
- [References](#References)


## Installation

This pipeline will require [nextflow](https://www.nextflow.io/docs/latest/getstarted.html). The rest fo the dependencies can be installed via [miniconda](https://github.com/isugifNF/polishCLR/tree/main#miniconda), [Docker](https://github.com/isugifNF/polishCLR/tree/main#docker), or [Singularity](https://github.com/isugifNF/polishCLR/tree/main#singularity)

```
nextflow -version

nextflow run isugifNF/polishCLR -r main \
  --help
```

For ag100pest projects, we have been installing dependencies in a [miniconda](https://docs.conda.io/en/latest/miniconda.html) environment. Dependencies may also be installed via [docker](https://www.docker.com/) or [singularity](https://sylabs.io/guides/3.0/user-guide/quick_start.html).

To update to the latest version of the pipeline, run:

```
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

<!-- # fill this in
If you are on MacOS with a Silicon chip (M1 or M2) you may run into the following warning.

```
WARNING: The requested image's platform (linux/amd64) does not match the detected host platform (linux/arm64/v8) and no specific platform was requested
```

You will need to add the following flag to specify the platform.
-->

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

## Basic Run

Print usage statement

  ```
  nextflow run isugifNF/polishCLR -r main --help

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
   --mitochondrial_assembly       mitochondrial assembly will be concatenated to the assemblies before polishing [default: false]

   Either FALCON (or FALCON Unzip) assembly:
   --primary_assembly             genome assembly fasta file to polish
   --alternate_assembly           if alternate/haplotig assembly file is provided, will be concatenated to the primary assembly before polishing [default: false]
   --falcon_unzip                 if primary assembly has already undergone falcon unzip [default: false]. If true, will Arrow polish once instead of twice.

   Or TrioCanu assembly
   --paternal_assembly            paternal genome assembly fasta file to polish
   --maternal_assembly            maternal genome assembly fasta file to polish

   Pick Step 1 (arrow, purgedups) or Step 2 (arrow, freebayes, freebayes)
   --step                         Run step 1 or step 2 (default: 1)

   Optional modifiers
   --species                      if a string is given, rename the final assembly by species name [default:false]
   --k                            kmer to use in MerquryQV scoring [default:21]
   --same_specimen                if illumina and pacbio reads are from the same specimen [default: true].
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

Example input data for each of the three input cases from _Helicoverpa zea_ on Ag Data Commons: https://data.nal.usda.gov/dataset/data-polishclr-example-input-genome-assemblies. SRAs from NCBI for pacbio and illumina data are available through the BioProject: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA804956

```
# If you have installed dependencies in a miniconda environment
source activate ${PWD}/env/polishCLR_env

SP_DIR=/project/ag100pest/Helicoverpa_zea_male

## === Step 1 (up to purge_dups)
nextflow run isugifNF/polishCLR -r main \
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
 nextflow run isugifNF/polishCLR -r main \
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
Key Outputs are found in
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

[timeline.html](https://isugifnf.github.io/polishCLR/timeline.html) | [report.html](https://isugifnf.github.io/polishCLR/report.html)

# References

Rhie, A., McCarthy, S. A., Fedrigo, O., Damas, J., Formenti, G., Koren, S., Uliano-Silva, M., Chow, W., Fungtammasan, A., Kim, J., Lee, C., Ko, B. J., Chaisson, M., Gedman, G. L., Cantin, L. J., Thibaud-Nissen, F., Haggerty, L., Bista, I., Smith, M., . . . Jarvis, E. D. (2021). Towards complete and error-free genome assemblies of all vertebrate species. Nature, 592(7856), 737-746. https://doi.org/10.1038/s41586-021-03451-0 

Rhie, A., Walenz, B. P., Koren, S., & Phillippy, A. M. (2020). Merqury: reference-free quality, completeness, and phasing assessment for genome assemblies. Genome Biol, 21(1), 245. https://doi.org/10.1186/s13059-020-02134-9.

Formenti, G., Rhie, A., Balacco, J., Haase, B., Mountcastle, J., Fedrigo, O., Brown, S., Capodiferro, M. R., Al-Ajli, F. O., Ambrosini, R., Houde, P., Koren, S., Oliver, K., Smith, M., Skelton, J., Betteridge, E., Dolucan, J., Corton, C., Bista, I., . . . Vertebrate Genomes Project, C. (2021). Complete vertebrate mitogenomes reveal widespread repeats and gene duplications. Genome Biol, 22(1), 120. https://doi.org/10.1186/s13059-021-02336-9 
