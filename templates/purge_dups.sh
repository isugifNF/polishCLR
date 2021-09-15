#! /usr/bin/env bash
# === Inputs
# genome_fasta=genome assembly file (cns_p_ctg.fasta)
# haplo_fasta=haplotype genome (cns_h_ctg.fasta)
# pacbio_reads=*subreads.fasta
# === Outputs

module load minimap2
module load python_3
#module load purge_dups # not sure how module differs from most recent git
export PATH="/project/ag100pest/software/purge_dups/bin/:$PATH"

PROC=\$((`nproc`))

${samtools_app} faidx ${primary_assembly}

${minimap2_app} -xmap-pb -t \${PROC}  ${primary_assembly} ${pacbio_reads} | \
  ${gzip_app} -c -  > p_mapping.paf.gz

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

${minimap2_app} -xmap-pb -t \${PROC}  h_${haplo_fasta} ${pacbio_reads} |\
  ${gzip_app} -c -  > h_mapping.paf.gz

${pbcstat_app} h_mapping.paf.gz

${calcuts_app} PB.stat > h_cutoffs 2> h_calcuts.log

${split_fa_app} h_${haplo_fasta} > h_${haplo_fasta}.split

${minimap2_app} -xasm5 -DP -t \${PROC} h_${haplo_fasta}.split h_${haplo_fasta}.split |\
  ${gzip_app} -c - > h_${haplo_fasta}.split.self.paf.gz

${purge_dups_app} -2 -T h_cufoffs -c PB.base.cov h_${haplo_fasta}.split.self.paf.gz > h_dups.bed 2> h_purge_dups.log

${get_seqs_app} -e h_dups.bed h_${haplo_fasta} -p haps

### TODO pull this out as a separate process in nextflow
echo "Purged alternate, running bbtools stats.sh on each assembly"
module load bbtools
## On Atlas
#export PATH="/project/ag100pest/software/bbmap/:$PATH"
stats.sh -Xmx2048m primary.purged.fa > primary_purged.fa.stats
stats.sh -Xmx2048m primary.hap.fa > primary_hap.fa.stats
stats.sh -Xmx2048m haps.purged.fa > haps_purged.fa.stats
stats.sh -Xmx2048m haps.hap.fa > haps_hap.fa.stats

## rename to play nicely with nextflow shortname 
mv primary.hap.fa primary_hap.fa
mv primary.purged.fa primary_purged.fa
mv haps.hap.fa haps_hap.fa
mv haps.purged.fa haps_purged.fa
#mv h_${haplo_fasta} h_genome.fa
