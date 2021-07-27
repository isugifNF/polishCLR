#! /usr/bin/env bash

# lookup_table.txt from bin folder, maybe need to run genomescope?
merfin -polish \
  -sequence ${genome_fasta} \
  -readmers ${meryldb} \
  -peak ${peak} \
  -prob lookup_table.txt \
  -vcf ${vcf} \
  -output ${i}_merfin
