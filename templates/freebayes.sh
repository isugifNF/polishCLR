#! /usr/bin/env bash
${freebayes_app} \
  --region "${window}" \
  --min-mapping-quality 0 \
  --min-coverage 3 \
  --min-supporting-allele-qsum 0 \
  --ploidy 2 \
  --min-alternate-fraction 0.2 \
  --max-complex-gap 0 \
  --bam ${illumina_bam} \
  --vcf ${illumina_bam.simpleName}_${window.replace(':','_').replace('|','_')}.vcf \
  --fasta-reference ${assembly_fasta}
