#! /usr/bin/env bash
# === Inputs
# genome_fasta=genome assembly file (cns_p_ctg.fasta)
# haplo_fasta=haplotype genome (cns_h_ctg.fasta)
# pacbio_reads=*subreads.fasta
# === Outputs

module load python_3

## These should all be included in the environment.yml
# module load minimap2
# module load samtools
# module load bamtools 
# module load purge_dups # not sure how module differs from most recent git

## Not sure if we need this still - testing
# export PATH="/project/ag100pest/software/purge_dups/bin/:$PATH"

PROC=\$((`nproc`))

${samtools_app} faidx ${primary_assembly}

for x in `${pacbio_reads} |  sed 's/\.bam//g'`; do
       ${minimap2_app} -xmap-pb $primary_assembly ${x}.fasta | $gzip_app -c - > p_mapping.paf.gz
done

${pbcstat_app} p_mapping.paf.gz
${hist_plot_py} PB.stat primary_hist

${calcuts_app} PB.stat > p_cutoffs 2> p_calcuts.log 

${split_fa_app} ${primary_assembly} > ${primary_assembly}.split

${minimap2_app} -xasm5 -DP -t \${PROC} ${primary_assembly}.split ${primary_assembly}.split | gzip -c - > ${primary_assembly}.split.self.paf.gz

${purge_dups_app} -2 -T p_cufoffs -c PB.base.cov ${primary_assembly}.split.self.paf.gz > p_dups.bed 2> p_purge_dups.log

## -e flag to only remove haplotypic regions at the ends of contigs
${get_seqs_app} -e p_dups.bed ${primary_assembly} -p primary

echo "Purged primary, onto the alternate"

## do the same with the "alternate" haps after cating
cat primary.hap.fa  ${haplo_fasta} >  h_${haplo_fasta}

${samtools_app} faidx h_${haplo_fasta}

for x in `${pacbio_reads} |  sed 's/\.bam//g'`; do
       ${minimap2_app} -xmap-pb h_${haplo_fasta} ${x}.fasta | $gzip_app -c - > h_mapping.paf.gz
done

${pbcstat_app} h_mapping.paf.gz

${hist_plot_py} PB.stat h_hist

${calcuts_app} PB.stat > h_cutoffs 2> h_calcuts.log

${split_fa_app} h_${haplo_fasta} > h_${haplo_fasta}.split

${minimap2_app} -xasm5 -DP -t \${PROC} h_${haplo_fasta}.split h_${haplo_fasta}.split |\
  ${gzip_app} -c - > h_${haplo_fasta}.split.self.paf.gz

${purge_dups_app} -2 -T h_cufoffs -c PB.base.cov h_${haplo_fasta}.split.self.paf.gz > h_dups.bed 2> h_purge_dups.log

${get_seqs_app} -e h_dups.bed h_${haplo_fasta} -p haps

echo "Purged alternate, running bbtools stats.sh on each assembly"

## On Atlas
#export PATH="/project/ag100pest/software/bbmap/:$PATH"
stats.sh -Xmx2048m primary.purged.fa > primary_purged.fa.stats
stats.sh -Xmx2048m primary.hap.fa > primary_hap.fa.stats
stats.sh -Xmx2048m haps.purged.fa > haps_purged.fa.stats
stats.sh -Xmx2048m haps.hap.fa > haps_hap.fa.stats

## rename to play nicely with nextflow shortname 
mv primary.purged.fa primary_purged.fa
mv haps.purged.fa haps_purged.fa
