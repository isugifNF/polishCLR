#! /usr/bin/env bash
PROC=\$(((`nproc`-1)*3/4+1))
${gcpp_app} --algorithm=arrow \
  -x 10 -X 120 -q 0 \
  -j \${PROC} -w "$window" \
  -r ${assembly_fasta} ${pacbio_bam} \
  -o ${window.replace(':','_').replace('|','_')}.vcf,${window.replace(':','_').replace('|','_')}.fasta
