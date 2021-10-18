#! /usr/bin/env bash
PROC=\$(((`nproc`-1)*3/4+1))
PROC2=\$(((`nproc`-1)*1/4+1))
mkdir tmp

# for multiple pacbio subread files
ls ${pacbio_read} > bam.fofn

${pbmm2_app} align -j \$PROC ${assembly_fasta} bam.fofn | \
  ${samtools_app} sort -T tmp -m 8G --threads \$PROC2 - > ${assembly_fasta.simpleName}_aln.bam
${samtools_app} index -@ \${PROC} ${assembly_fasta.simpleName}_aln.bam
