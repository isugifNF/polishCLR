#! /usr/bin/env bash
# === Inputs
# genome_fasta=genome assembly file (cns_p_ctg.fasta)
# haplo_fasta=haplotype genome (cns_h_ctg.fasta)
# pacbio_reads=*subreads.fasta
# === Outputs


PROC=\$((`nproc`))

${samtools_app} faidx ${primary_assembly}

for x in ${pacbio_reads.simpleName} ; do
       ${minimap2_app} -xmap-pb $primary_assembly \${x}.fasta | $gzip_app -c - > \${x}_p_mapping.paf.gz
done

${pbcstat_app} *_p_mapping.paf.gz

${calcuts_app} PB.stat > p_cutoffs 2> p_calcuts.log 

#${hist_plot_py} PB.stat --cutoffs p_cutoffs primary_hist

${split_fa_app} ${primary_assembly} > ${primary_assembly}.split

${minimap2_app} -xasm5 -DP -t \${PROC} ${primary_assembly}.split ${primary_assembly}.split | gzip -c - > ${primary_assembly}.split.self.paf.gz

${purge_dups_app} -2 -T p_cufoffs -c PB.base.cov ${primary_assembly}.split.self.paf.gz > p_dups.bed 2> p_purge_dups.log

## -e flag to only remove haplotypic regions at the ends of contigs
${get_seqs_app} -e p_dups.bed ${primary_assembly} -p primary

echo "Purged primary, onto the alternate"

## do the same with the "alternate" haps after cating
cat primary.hap.fa  ${haplo_fasta} >  h_${haplo_fasta}

${samtools_app} faidx h_${haplo_fasta}

for x in ${pacbio_reads.simpleName} ; do
       ${minimap2_app} -xmap-pb $primary_assembly \${x}.fasta | $gzip_app -c - > \${x}_h_mapping.paf.gz
done

${pbcstat_app} *_h_mapping.paf.gz

${calcuts_app} PB.stat > h_cutoffs 2> h_calcuts.log

#${hist_plot_py} PB.stat --cutoffs h_cutoffs h_hist

${split_fa_app} h_${haplo_fasta} > h_${haplo_fasta}.split

${minimap2_app} -xasm5 -DP -t \${PROC} h_${haplo_fasta}.split h_${haplo_fasta}.split |\
  ${gzip_app} -c - > h_${haplo_fasta}.split.self.paf.gz

${purge_dups_app} -2 -T h_cufoffs -c PB.base.cov h_${haplo_fasta}.split.self.paf.gz > h_dups.bed 2> h_purge_dups.log

${get_seqs_app} -e h_dups.bed h_${haplo_fasta} -p haps

## create and add kat module 

## rename to play nicely with nextflow simplename 
mv primary.purged.fa primary_purged.fa
mv haps.purged.fa haps_purged.fa
