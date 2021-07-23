#! /usr/bin/env bash
PROC=\$(((`nproc`-1)/2+1))
${parallel_app} -j 2 "${bzcat_app} {1} | ${pigz_app} -p \${PROC} > {1/.}.gz" ::: *.bz2
