#! /usr/bin/env bash
cat ${vcfs.get(0)} | grep "^#" > ${i}_consensus.vcf
cat ${vcfs} | grep -v "^#" >> ${i}_consensus.vcf
