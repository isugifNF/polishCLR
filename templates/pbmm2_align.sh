#! /usr/bin/env bash
PROC=\$(((`nproc`-1)*3/4+1))
PROC2=\$(((`nproc`-1)*1/4+1))
mkdir tmp
${pbmm2_app} align -j \$PROC ${assembly_fasta} ${pacbio_read} | \
  ${samtools_app} sort -T tmp -m 8G --threads \$PROC2 - > ${pacbio_read.simpleName}_aln.bam
${samtools_app} index -@ \${PROC} ${pacbio_read.simpleName}_aln.bam
