#! /usr/local/bin bash

${samtools_app} faidx ${genome_fasta}
grep ">pri_" ${genome_fasta} | cut -f1 > pri.list
${samtools_app} faidx -r pri.list ${genome_fasta} > p_${genome_fasta}

grep ">mit_" ${genome_fasta} | cut -f1 > mit.list
${samtools_app} faidx -r mit.list ${genome_fasta} > m_${genome_fasta}

grep ">alt_" ${genome_fasta} | cut -f1 > alt.list
${samtools_app} faidx -r alt.list ${genome_fasta} > a_${genome_fasta}