#! /usr/bin/env bash
PROC=\$((`nproc`))

# == reshape_arrow.sh
# https://github.com/arangrhie/merfin/blob/master/scripts/reformat_arrow/reshape_arrow.sh
output=${vcf.baseName}. # assumes input.vcf
grep -v "#" \${output}.vcf | sed 's/,/;/g' > \${output}.temp.reshaped.vcf
${bcftools_app} view -h \${output}.vcf > \${output}.temp.reshaped.header.vcf
cat \${output}.temp.reshaped.header.vcf \${output}.temp.reshaped.vcf > \${output}.temp.reshaped.combined.vcf
rm \${output}.temp.reshaped.header.vcf \${output}.temp.reshaped.vcf
${bcftools_app} annotate -h extra_header.vcf \${output}.temp.reshaped.combined.vcf > \${output}.temp.reshaped.vcf
${bcftools_app} view -h \${output}.temp.reshaped.vcf | sed 's/\tINFO/\tINFO\tFORMAT\tIND/g' > \${output}.reshaped.vcf
rm \${output}.temp.reshaped.vcf 
${bcftools_app} view -H \${output}.temp.reshaped.combined.vcf | awk -F"\t" -v OFS="\t" '{gsub(/DP=/,".\tGT:DP\t1/1:",\$8);print \$0}' >> \${output}.reshaped.vcf
# https://github.com/arangrhie/merfin/wiki/Best-practices-for-Merfin#preparing-input-vcf
bcftools view -Oz -i ‘(GT="AA" || GT="Aa")' \${output}.reshaped.vcf > out.vcf.gz
${bcftools_app} view out.vcf.gz -Oz > \${output}.reshaped.vcf.gz
rm \${output}.reshaped.vcf 
rm \${output}.temp.reshaped.combined.vcf
rm out.vcf
