#! /usr/bin/env bash

# Calculate histogram
${meryl_app} histogram ${illumina_meryl} > illumina.hist
awk '\$1>5 {print}' illumina.hist | sort -k 2n | tail -n 1 | awk '{print \$1}' > peak.txt
