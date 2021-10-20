#! /usr/bin/env bash
${bcftools_app} view -Oz ${vcf} > ${vcf}.gz
${bcftools_app} index ${vcf}.gz
${bcftools_app} consensus ${vcf}.gz -f ${genome_fasta} -Hla > ${outdir}_consensus.fasta
