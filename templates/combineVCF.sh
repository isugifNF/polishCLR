#! /usr/bin/env bash

OUTVCF=${i}_consensus.vcf

# Merge by sections (1) headers up to contig, (2) all contig headers, (3) headers post contigs, (4) snp data
cat ${vcfs.get(0)} | sed -n '1,/##reference=/'p > \$OUTVCF
cat ${vcfs} | grep -h "##contig=" >> \$OUTVCF
cat ${vcfs.get(0)} | sed -n '/##INFO=/,/#CHROM/'p >> \$OUTVCF
cat ${vcfs} | grep -hv "#" >> \$OUTVCF

# === The grep statements only pull one contig header
#cat ${vcfs.get(0)} | grep "^#" > ${i}_consensus.vcf
#cat ${vcfs} | grep -v "^#" >> ${i}_consensus.vcf
# === Was going to use the following instead, but hit a duplicate sample name error
#PROC=\$((`nproc`))
#${parallel_app} -j \$PROC "bcftools view -Oz {1} > {1}.gz" ::: ${vcfs}
#${parallel_app} -j \$PROC "bcftools index {1}.gz" ::: ${vcfs}
#${bcftools_app} merge --force-samples *vcf.gz > ${i}_consensus.vcf
