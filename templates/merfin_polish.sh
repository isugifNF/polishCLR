#! /usr/bin/env bash

OUTNAME=`echo $outdir | sed 's:/:_:g'`

${merfin_app} -polish \
  -sequence ${genome_fasta} \
  -seqmers ${genome_meryl} \
  -readmers ${meryldb} \
  -peak ${peak} \
  -vcf ${vcf} \
  -output \${OUTNAME}_merfin
