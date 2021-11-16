#! /usr/bin/env bash

OUTNAME=`echo $outdir | sed 's:/:_:g'`

${bcftools_app} view -Oz ${vcf} > ${vcf}.gz
${bcftools_app} index ${vcf}.gz
${bcftools_app} consensus ${vcf}.gz -f ${genome_fasta} -Hla > \${OUTNAME}_consensus.fasta
