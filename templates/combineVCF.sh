#! /usr/bin/env bash
PROC=\$((`nproc`))
${parallel_app} -j \$PROC "bcftools view -Oz {1} > {1}.gz" ::: ${vcfs}
${parallel_app} -j \$PROC "bcftools index {1}.gz" ::: ${vcfs}
${bcftools_app} merge *vcf.gz > ${i}_consensus.vcf