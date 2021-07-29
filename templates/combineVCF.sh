#! /usr/bin/env bash

cat ${vcfs.get(0)} | grep "^#" > ${i}_consensus.vcf
cat ${vcfs} | grep -v "^#" >> ${i}_consensus.vcf

# === Was going to use the following instead, but hit a duplicate sample name error
#PROC=\$((`nproc`))
#${parallel_app} -j \$PROC "bcftools view -Oz {1} > {1}.gz" ::: ${vcfs}
#${parallel_app} -j \$PROC "bcftools index {1}.gz" ::: ${vcfs}
#${bcftools_app} merge --force-samples *vcf.gz > ${i}_consensus.vcf
