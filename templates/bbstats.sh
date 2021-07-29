#! /usr/bin/env bash
# Desc: Prints out length distibutions, GC, etc of each assembly, could be added to pipeline at the end

module load bbtools

echo "Assmbly stats of $assembly_fasta  according to bbtools stats.sh"

stats.sh in=$assembly_fasta out=${assembly_fasta.simpleName}
