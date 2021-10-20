#! /usr/bin/env bash

${merfin_app} -polish \
  -sequence ${genome_fasta} \
  -seqmers ${genome_meryl} \
  -readmers ${meryldb} \
  -peak ${peak} \
  -vcf ${vcf} \
  -output ${outdir}_merfin
