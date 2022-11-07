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
nextflow run isugifNF/polishCLR -r main  \
  --primary_assembly "data/primary.fasta" \
  --mitocondrial_assembly "data/mitochondrial.fasta" \
  --illiumina_reads "data/illumina/*_{R1,R2}.fasta.bz" \
  --pacbio_reads "data/pacbio/pacbio.subreads.bam" \
  --step 1 \
  --falcon-unzip false \
  -profile slurm
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
