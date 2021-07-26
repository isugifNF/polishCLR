#! /usr/bin/env bash
${bcftools_app} view -Oz ${vcf} > ${vcf.simpleName}.gz
${bcftools_app} index ${vcf.simpleName}.gz
${bcftools_app} consensus ${vcf.simpleName}.gz -f ${genome_fasta} -H 1 > consensus_04.fasta
