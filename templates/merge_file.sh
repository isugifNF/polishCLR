#! /usr/local/bin bash
cat ${asm_file} | sed 's/>/>pri_/g' > ${asm_file.simpleName}_mito.fasta
echo "" >> ${asm_file.simpleName}_mito.fasta
cat ${mito_file} | sed 's/>/>mit_/g' >> ${asm_file.simpleName}_mito.fasta
echo "" >> ${asm_file.simpleName}_mito.fasta
cat ${alt_file} | sed 's/>/>alt_/g' >> ${asm_file.simpleName}_mito.fasta
cat ${asm_file.simpleName}_mito.fasta | grep -v "^\$" > temp.txt
mv temp.txt ${asm_file.simpleName}_mito.fasta