#! /usr/bin/env bash

printf "======================================================="
printf "uniq kmers in asm | kmers in both asm and reads | QV | Error rate"
printf "======================================================="
${merqury_sh} $illumina_db $assembly_fasta ${assembly_fasta.simpleName}
