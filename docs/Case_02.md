---
title: "Case 2: Haplotype-Resolved"
sidebar: toc
maxdepth: 1
sort: 2
---

# Case 2: Haplotype-Resolved

This describes Case 2 of the three input cases:

* Case 1: An unresolved primary assembly with associated contigs (the output of FALCON 2-asm) or without (e.g., the output of Canu or wtdbg2).
* **Case 2: A haplotype-resolved but unpolished set (e.g., the output of FALCON-Unzip 3-unzip).**
* Case 3: IDEAL! A haplotype-resolved, CLR long-read, Arrow-polished set of primary and alternate contigs (e.g., the output of FALCON-Unzip 4-polish).

### Recommended parameters

```
nextflow run isugifNF/polishCLR --main \
  --primary_assembly "data/primary.fasta" \
  --alternate_assembly "data/alternate.fasta" \
  --mitocondrial_assembly "data/mitochondrial.fasta" \
  --illumina_reads "data/illumina/*_{R1,R2}.fasta.bz" \
  --pacbio_reads "data/pacbio/pacbio.subreads.bam" \
  --step 1 \
  --arrow01 \
  -profile slurm \
  -resume
```
