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

${split_fa_app} ${primary_assembly} > ${primary_assembly}.split

${minimap2_app} -xasm5 -DP -t \${PROC} ${primary_assembly}.split ${primary_assembly}.split | gzip -c - > ${primary_assembly}.split.self.paf.gz

${purge_dups_app} -2 -T p_cufoffs -c PB.base.cov ${primary_assembly}.split.self.paf.gz > p_dups.bed 2> p_purge_dups.log

## -e flag to only remove haplotypic regions at the ends of contigs
${get_seqs_app} -e p_dups.bed ${primary_assembly} -p primary

echo "Purged primary, onto the alternate"

## do the same with the "alternate" haps after cating
${samtools_app} faidx primary.hap.fa

for x in ${pacbio_reads.simpleName} ; do
       ${minimap2_app} -xmap-pb  primary.hap.fa \${x}.fasta | $gzip_app -c - > \${x}_h_mapping.paf.gz
done

${pbcstat_app} *_h_mapping.paf.gz

${calcuts_app} PB.stat > h_cutoffs 2> h_calcuts.log

${split_fa_app} primary.hap.fa > primary.hap.split

${minimap2_app} -xasm5 -DP -t \${PROC} primary.hap.split primary.hap.split |\
  ${gzip_app} -c - > primary.hap.split.self.paf.gz

${purge_dups_app} -2 -T h_cufoffs -c PB.base.cov primary.hap.split.self.paf.gz > h_dups.bed 2> h_purge_dups.log

${get_seqs_app} -e h_dups.bed primary.hap.fa -p haps

## create and add kat module 

## rename to play nicely with nextflow simplename 
mv primary.purged.fa primary_purged.fa
mv haps.purged.fa haps_purged.fa
