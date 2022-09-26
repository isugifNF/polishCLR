#! /usr/bin/env bash
# === Inputs
# genome_fasta=genome assembly file from maternal or paternal haplotyopes
# pacbio_reads=*subreads.fasta from F1
# === Outputs
# primary_hist
# fastas

module load minimap2
module load python_3
export PATH="/project/ag100pest/software/purge_dups/bin/:$PATH"

PROC=\$((`nproc`))

## TODO make sure that the short names inlcude maternal/paternal 

${samtools_app} faidx ${primary_assembly}

for x in ${pacbio_reads.simpleName} ; do
       ${minimap2_app} -xmap-pb ${params.minimap2_params} $primary_assembly \${x}.fasta | $gzip_app -c - > \${x}_p_mapping.paf.gz
done

${pbcstat_app} *_${primary_assembly.shortName}_mapping.paf.gz
${hist_plot_py} PB.stat ${primary_assembly.shortName}_hist

${calcuts_app} PB.stat > p_cutoffs 2> p_calcuts.log 

${split_fa_app} ${primary_assembly} > ${primary_assembly.shortName}.split

${minimap2_app} -xasm5 ${params.minimap2_params} -DP -t \${PROC} ${primary_assembly.shortName}.split ${primary_assembly.shortName}.split | gzip -c - > ${primary_assembly.shortName}.split.self.paf.gz

${purge_dups_app} ${purge_dups_params} -c PB.base.cov ${primary_assembly.shortName}.split.self.paf.gz > ${primary_assembly.shortName}_dups.bed 2> ${primary_assembly.shortName}_purge_dups.log

## In trio assemblies there shouldn't be haplotigs, so remove these from the bed file output by purge_dups.bed
grep 'JUNK\|OVLP' ${primary_assembly.shortName}_dups.bed > ${primary_assembly.shortName}_dups_JUNK_OVLP.bed

## -e flag to only remove haplotypic regions at the ends of contigs
${get_seqs_app} -e dups_JUNK_OVLP.bed ${primary_assembly} -p ${primary_assembly.shortName}

module load bbtools
stats.sh -Xmx2048m ${primary_assembly.shortName}.purged.fa > ${primary_assembly.shortName}_purged.fa.stats
stats.sh -Xmx2048m ${primary_assembly.shortName}.hap.fa > ${primary_assembly.shortName}_hap.fa.stats

## rename to play nicely with nextflow shortname 
mv ${primary_assembly.shortName}.hap.fa ${primary_assembly.shortName}_hap.fa
mv ${primary_assembly.shortName}.purged.fa ${primary_assembly.shortName}_purged.fa

