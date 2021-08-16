#! /usr/bin/env bash

# lookup_table.txt from bin folder, maybe need to run genomescope?
${merfin_app} -polish \
  -sequence ${genome_fasta} \
  -seqmers ${genome_meryl} \
  -readmers ${meryldb} \
  -peak ${peak} \
  -vcf ${vcf} \
  -output ${i}_merfin
# merfin -polish \
#   -sequence GENOME_FASTA \
#   -seqmers genome.meryl \
#   -readmers ILLUMINA_MERYL \
#   -peak 79 \
#   -vcf arrow_merged.reshaped.vcf.gz \
#   -output arrow_merfinpolish.gz > arrow_merfinpolish.out

