#! /usr/bin/env bash
PROC=\$((`nproc))

# == merge vcf
${parallel_app} -j \$PROC "bcftools view -Oz {1} > {1}.gz" ::: ${vcf}
${parallel_app} -j \$PROC "bcftools index {1}" ::: ${vcf}
bcftools merge *.gz > ${i}_merged.vcf
grep "#" ${i}_merged.vcf > ${i}_reshaped.vcf
grep -v "#" ${i}_merged.vcf | sed 's/,/;/g' >> ${i}_reshaped.vcf

# -- preprocess arrow
#bcftools view -h ${vcf} > ${vcf.simpleName}.temp.reshaped.header.vcf
#cat ${vcf.simpleName}.temp.reshaped.header.vcf ${vcf.simpleName}.temp.reshaped.vcf > ${vcf.simpleName}.temp.reshaped.combined.vcf
#rm ${vcf.simpleName}.temp.reshaped.header.vcf ${vcf.simpleName}.temp.reshaped.vcf
#bcftools annotate -h extra_header.vcf ${vcf.simpleName}.temp.reshaped.combined.vcf > ${vcf.simpleName}.temp.reshaped.vcf
#bcftools view -h ${vcf.simpleName}.temp.reshaped.vcf | sed 's/\tINFO/\tINFO\tFORMAT\tIND/g' > ${vcf.simpleName}.reshaped.vcf
#rm ${vcf.simpleName}.temp.reshaped.vcf 
#bcftools view -H ${vcf.simpleName}.temp.reshaped.combined.vcf | awk -F"\\t" -v OFS="\\t" '{gsub(/DP=/,".\\tGT:DP\\t1/1:",\$8);print \$0}' >> ${vcf.simpleName}.reshaped.vcf
#bcftools view ${vcf.simpleName}.reshaped.vcf -Oz > ${vcf.simpleName}.reshaped.vcf.gz
#rm ${vcf.simpleName}.reshaped.vcf 
#rm ${vcf.simpleName}.temp.reshaped.combined.vcf
