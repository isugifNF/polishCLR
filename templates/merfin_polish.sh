#! /usr/bin/env bash




${merfin_app} -polish \
  -sequence ${genome_fasta} \
  -readmers ${meryldb} \
  -peak ${peak} \
  -prob lookup_table.txt \
  -vcf ${vcf} \
  -output ${i}_merfin
