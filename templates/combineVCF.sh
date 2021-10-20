#! /usr/bin/env bash

OUTVCF=${outdir}_consensus.vcf

# Merge by sections (1) headers up to contig, (2) all contig headers, (3) headers post contigs, (4) snp data
cat ${vcfs.get(0)} | sed -n '1,/##reference=/'p > \$OUTVCF
cat ${vcfs} | grep -h "##contig=" >> \$OUTVCF
cat ${vcfs.get(0)} | sed -n '/##INFO=/,/#CHROM/'p >> \$OUTVCF
cat ${vcfs} | grep -hv "#" >> \$OUTVCF

