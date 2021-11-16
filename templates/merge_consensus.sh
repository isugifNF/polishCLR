#! /usr/bin/env bash

OUTNAME=`echo "$outdir" | sed 's:/:_:g'`
cat ${windows_fasta} > \${OUTNAME}_consensus.fasta
