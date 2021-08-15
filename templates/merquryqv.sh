#! /usr/bin/env bash

printf "======================================================="
printf "uniq kmers in asm | kmers in both asm and reads | QV | Error rate"
printf "======================================================="
${merqury_sh} $illumina_db $assembly_fasta ${assembly_fasta.simpleName}
printf "=======================================================" > ${assembly_fasta.simpleName}_qv.txt
printf "uniq kmers in asm | kmers in both asm and reads | QV | Error rate" >> ${assembly_fasta.simpleName}_qv.txt
printf "=======================================================" >> ${assembly_fasta.simpleName}_qv.txt
cat  ${assembly_fasta.simpleName}.qv >> ${assembly_fasta.simpleName}_qv.txt

# == Get single QV value
cat ${assembly_fasta.simpleName}.qv | awk -F'\t' '{print \$4}' > merqury.qv