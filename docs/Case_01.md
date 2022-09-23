---
title: "Case 1: Primary Assembly"
sidebar: toc
maxdepth: 1
sort: 1
---

# Case 1
 **A primary assembly without haplotype  (e.g., the output of Canu or wtdbg2).**

Step 1 runs another round of Arrow polishing and then purge_dups to split off alternate haplotypes.

```
nextflow run <path/to/polishCLR>/main.nf  \
  --primary_assembly "data/primary.fasta" \
  --mitocondrial_assembly "data/mitochondrial.fasta" \
  --illiumina_reads "data/illumina/*_{R1,R2}.fasta.bz" \
  --pacbio_reads "data/pacbio/pacbio.subreads.bam" \
  --step 1 \
  --falcon-unzip false \
  -profile slurm
```

Step 2 runs another round of Arrow polishing with the PacBio reads, then polishes with short-reads with two rounds of FreeBayes.

Provide the purged primary and alternate contigs from purge dups, and mitochondrial genome as input - or, if scaffolding data, like Hi-C, are available to you, you should scaffold the output of Step 1. Don't forget to include parameter flags `--step 2` and `resume` to this command. 

`haps_purged.fa` and `primary_purged.fa`

```
 nextflow run <path/to/polishCLR>/main.nf \
  --primary_assembly "primary_purged.fa" \
  --alternate_assembly "haps_purged.fa" \
  --mitochondrial_assembly "Otur_mtDNA_contig.fasta" \
  --illumina_reads "../RawPolishingData/JAMW*{R1,R2}.fastq.bz2" \
  --pacbio_reads "../RawSequelData/m*.subreads.bam" \
  --step 2 \
  -resume
  ```
