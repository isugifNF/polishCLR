#! /usr/bin/env bash
PROC=\$((`nproc`))

# == Jellyfish to get peak
${parallel_app} -j \$PROC "gunzip {1}" ::: $illumina_reads
jellyfish count -C -m "${params.k}" -t \$PROC -s 3000000000 *.fastq -o illumina.jf
jellyfish histo -t \$PROC illumina.jf > illumina.hist
cat illumina.hist | sort -k 2n | tail -n 1 |awk '{print \$1}' > peak.txt
