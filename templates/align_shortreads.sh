#! /usr/bin/env bash
PROC=\$(((`nproc`-1)*3/4+1))
PROC2=\$(((`nproc`-1)*1/4+1))
mkdir tmp
${bwamem2_app} index ${assembly_fasta}
${bwamem2_app} mem -SP -t \$PROC ${assembly_fasta} ${illumina_one} ${illumina_two} | \
    ${samtools_app} sort -T tmp -m 8G --threads \$PROC2 - > ${illumina_one.simpleName}_aln.bam
${samtools_app} index -@ \${PROC} ${illumina_one.simpleName}_aln.bam
