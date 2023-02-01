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

**Example Dataset**

We have listed some example files to test the pipeline based on Chromosome 30 Hzea at https://data.nal.usda.gov/dataset/data-polishclr-example-input-genome-assemblies.

Case 2 will take primary assembly from the `FALCON/3-unzip` folder.

| Param | Files | Download link|
|:--|:--|:--
| `--primary_assembly` | "all_p_ctg.fasta" | [all_p_ctg.fasta](https://data.nal.usda.gov/system/files/all_p_ctg.fasta)|
| `--alternate_assembly` | "all_h_ctg.fasta" |[all_h_ctg.fasta](https://data.nal.usda.gov/system/files/all_h_ctg.fasta)|
| `--mitochondrial_assembly` | "GCF_022581195.2_ilHelZeax1.1_mito.fa" | [GenBank download fasta](https://www.ncbi.nlm.nih.gov/nuccore/NC_061507.1?report=fasta)|
| `--illumina_reads` |"testpolish_{R1,R2}.fastq" | [testpolish_R1.fastq](https://data.nal.usda.gov/system/files/testpolish_R1.fastq), [testpolish_R2.fastq](https://data.nal.usda.gov/system/files/testpolish_R2.fastq) |
| `--pacbio_reads` | "test.1.filtered.bam" | [test.1.filtered.bam_.gz](https://data.nal.usda.gov/system/files/test.1.filtered.bam_.gz)|

### Recommended parameters

```
nextflow run isugifNF/polishCLR -r main \
  --primary_assembly "all_p_ctg.fasta" \
  --alternate_assembly "all_h_ctg.fasta" \
  --mitochondrial_assembly "GCF_022581195.2_ilHelZeax1.1_mito.fa" \
  --illumina_reads "*_{R1,R2}.fastq" \
  --pacbio_reads "test.1.filtered.bam_.gz" \
  --step 1 \
  --arrow01 \
  -profile slurm \
  -resume
```

**Note:** On some browsers, the dashes (-) and underscores (_) can be copied incorrectly.  So if you run into an error that says `not valid in the pipeline` try manually retyping those parameters.

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
