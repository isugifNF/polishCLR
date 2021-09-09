#! /usr/bin/env bash
# genome_fasta=genome assembly file (cns_p_ctg.fasta)
# haplotype_fasta=haplotype genome (cns_h_ctg.fasta)
# pacbio_reads=*subreads.fasta

PROC=\$((`nproc`))

${samtools_app} faidx ${genome_fasta}
${minimap2_app} -xmap-pb -t \${PROC} ${genome_fasta} ${pacbio_reads} |\
  ${gzip_app} -c - > p_mapping.paf.gz

${pbcstat_app} p_mapping.paf.gz

${hist_plot_py} PB.stat primary_hist

${calcuts_app} PB.stat > p_cutoffs 2> p_calcuts.log
${split_fa} p_genome.fasta > p_genome.fasta.split

${minimap2_app} -xasm5 -DP -t \$PROC p_genome.fasta.split p_genome.fasta.split |\
  ${gzip_app} -c - > p_genome.fasta.self.paf.gz

${purge_dups_app} -2 -T p_cufoffs -c PB.base.cov p_genome.fasta.split.self.paf.gz > p_dups.bed 2> p_purge_dups.log 

## -e flag to only remove haplotypic regions at the ends of contigs (repeat of above)
${get_seqs_app} -e p_dups.bed p_genome.fasta -p primary

## do the same with the "alternate" haps after cating 
${samtools_app} faidx ${haplotype_fasta}

${minimap2_app} -xmap-pb -t 40 ${haplotype_fasta} ${pacbio_reads} |\
 ${gzip_app} -c -  > h_mapping.paf.gz

${pbcstat_app} h_mapping.paf.gz

${calcuts_app} PB.stat > h_cutoffs 2> h_calcuts.log

${split_fa_app} ${haplotype_fasta} > h_genome.fasta.split

${minimap2_app} -xasm5 -DP -t 40 h_genome.fasta.split h_genome.fasta.split | gzip -c - > h_genome.fasta.split.self.paf.gz

${purge_dups_app} -2 -T h_cufoffs -c PB.base.cov h_genome.fasta.split.self.paf.gz > h_dups.bed 2> h_purge_dups.log

${get_seqs_app} -e h_dups.bed h_genome.fasta -p haps

# after done, merge both

stats.sh -Xmx2048m primary.purged.fa > primary.purged.fa.stats

stats.sh -Xmx2048m primary.hap.fa > primary.hap.fa.stats

stats.sh -Xmx2048m haps.purged.fa > haps.purged.fa.stats

stats.sh -Xmx2048m haps.hap.fa > haps.hap.fa.stats


