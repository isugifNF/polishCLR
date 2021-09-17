#! /usr/bin/env bash

# === Inputs
# genome_fasta = primary_assembly_merged.fasta
# === Outputs
# p_${genome_fasta}     # primary assembly
# a_${genome_fasta}     # alternative assembly
# m_${genome_fasta}     # mitochondrial assembly

${samtools_app} faidx ${genome_fasta}
grep ">pri_" ${genome_fasta} | cut -f1 | sed 's/>//g' > pri.list
${samtools_app} faidx -r pri.list ${genome_fasta} > pat_${genome_fasta}

grep ">mit_" ${genome_fasta} | cut -f1 | sed 's/>//g' > mit.list
${samtools_app} faidx -r mit.list ${genome_fasta} > mit_${genome_fasta}