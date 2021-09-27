#! /usr/bin/env bash

# === Inputs
# primary_assembly = p_ctg.fasta     # From CANU, similar to Falcon but without alternative assembly
# mito_assembly = mt.fasta           # From vgpMito pipeline

# === Outputs
# ${primary_assembly.simpleName}_merged.fasta

cat ${primary_assembly} | sed 's/>/>pri_/g' > ${primary_assembly.simpleName}_temp.fasta
echo "" >> ${primary_assembly.simpleName}_temp.fasta
cat ${mito_assembly} | sed 's/>/>mit_/g' >> ${primary_assembly.simpleName}_temp.fasta
cat ${primary_assembly.simpleName}_temp.fasta | grep -v "^\$" > ${primary_assembly.simpleName}_merged.fasta
