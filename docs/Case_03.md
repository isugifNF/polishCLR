---
title: "Case 3: Polished"
sidebar: toc
maxdepth: 1
sort: 3
---

# Case 2: Polished

This describes Case 3 of the three input cases:

* Case 1: An unresolved primary assembly with associated contigs (the output of FALCON 2-asm) or without (e.g., the output of Canu or wtdbg2).
* Case 2: A haplotype-resolved but unpolished set (e.g., the output of FALCON-Unzip 3-unzip).
* **Case 3: IDEAL! A haplotype-resolved, CLR long-read, Arrow-polished set of primary and alternate contigs (e.g., the output of FALCON-Unzip 4-polish).**

### Recommended parameters

```
nextflow run isugifNF/polishCLR -r main \
  --primary_assembly "data/primary.fasta" \
  --alternate_assembly "data/alternate.fasta" \
  --mitocondrial_assembly "data/mitochondrial.fasta" \
  --illumina_reads "data/illumina/*_{R1,R2}.fasta.bz" \
  --pacbio_reads "data/pacbio/pacbio.subreads.bam" \
  --step 1 \
  -profile slurm \
  -resume
```

Step 2 runs another round of Arrow polishing with the PacBio reads, then polishes with short-reads with two rounds of FreeBayes. We broke these two steps into seperate phases to allow for manual scaffolding.

Provide the purged primary `primary_purged.fa` and alternate contigs `haps_purged.fa` from purge_dups, and mitochondrial genome `mitochondrial.fasta` as input to step 2. 

If scaffolding data, like Hi-C, are available to you, you should scaffold the `primary_purged.fa` and provide that as input for the  `--primary_assembly`. 

Regardless don't forget to include parameter flags `--step 2` and `resume` to this command. 

```
 nextflow run isugifNF/polishCLR -r main \
  --primary_assembly "primary_purged.fa" \
  --alternate_assembly "haps_purged.fa" \
  --mitochondrial_assembly "data/mitochondrial.fasta" \
  --illumina_reads "data/illumina/*_{R1,R2}.fasta.bz" \
  --pacbio_reads "../RawSequelData/m*.subreads.bam" \
  --step 2 \
  -profile slurm \
  -resume
  ```
