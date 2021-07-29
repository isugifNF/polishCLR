#! /usr/bin/env bash
${samtools_app} faidx ${assembly_fasta}
cat ${assembly_fasta}.fai | awk '{print \$1 ":0-" \$2}' > win.txt
