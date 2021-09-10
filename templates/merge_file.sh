#! /usr/local/bin bash

# === Inputs
# primary_assembly = p_ctg.fasta     # From FALCON or FALCON-Unzip 
# alternate_assembly = a_ctg.fasta
# mito_assembly = mt.fasta           # From vgpMito pipeline

# === Outputs
# ${primary_assembly.simpleName}_merged.fasta

cat ${primary_assembly} | sed 's/>/>pri_/g' > ${primary_assembly.simpleName}_temp.fasta
echo "" >> ${primary_assembly.simpleName}_temp.fasta
cat ${mito_assembly} | sed 's/>/>mit_/g' >> ${primary_assembly.simpleName}_temp.fasta
echo "" >> ${primary_assembly.simpleName}_temp.fasta
cat ${alternate_assembly} | sed 's/>/>alt_/g' >> ${primary_assembly.simpleName}_temp.fasta
cat ${primary_assembly.simpleName}_temp.fasta | grep -v "^\$" > ${primary_assembly.simpleName}_merged.fasta