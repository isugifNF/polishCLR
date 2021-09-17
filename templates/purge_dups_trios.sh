#! /usr/bin/env bash
# === Inputs
# genome_fasta=genome assembly file from mom or dad
# pacbio_reads=*subreads.fasta from mom or dad respectively
# === Outputs
# primary_hist
# fastas

module load minimap2
module load python_3
export PATH="/project/ag100pest/software/purge_dups/bin/:$PATH"

PROC=\$((`nproc`))

## TODO make sure that the short names inlcude maternal/paternal 

${samtools_app} faidx ${primary_assembly}

${minimap2_app} -xmap-pb -t \${PROC}  ${primary_assembly} ${pacbio_reads} | \
  ${gzip_app} -c -  > p_mapping.paf.gz

${pbcstat_app} p_mapping.paf.gz
${hist_plot_py} PB.stat primary_hist

${calcuts_app} PB.stat > p_cutoffs 2> p_calcuts.log 

${split_fa_app} ${primary_assembly} > ${primary_assembly}.split

${minimap2_app} -xasm5 -DP -t \${PROC} ${primary_assembly}.split ${primary_assembly}.split | gzip -c - > ${primary_assembly}.split.self.paf.gz

${purge_dups_app} -2 -T p_cufoffs -c PB.base.cov ${primary_assembly}.split.self.paf.gz > p_dups.bed 2> p_purge_dups.log

## In trio assemblies there shouldn't be haplotigs, so remove these from the bed file output by purge_dups.bed
grep 'JUNK\|OVLP' p_dups.bed > dups_JUNK_OVLP.bed

## -e flag to only remove haplotypic regions at the ends of contigs
${get_seqs_app} -e dups_JUNK_OVLP.bed ${primary_assembly} -p primary

### TODO pull this out as a separate process in nextflow
echo "Purged alternate, running bbtools stats.sh on each assembly"
module load bbtools
stats.sh -Xmx2048m primary.purged.fa > ${primary_assembly.simpleName}_primary_purged.fa.stats
stats.sh -Xmx2048m primary.hap.fa > ${primary_assembly.simpleName}_primary_hap.fa.stats

## rename to play nicely with nextflow shortname 
mv primary.hap.fa ${primary_assembly.simpleName}_primary_hap.fa
mv primary.purged.fa ${primary_assembly.simpleName}_primary_purged.fa

