#! /usr/bin/env nextflow

nextflow.enable.dsl=2

process PURGE_DUPS {
  publishDir "${params.outdir}/$outdir", mode:'copy'
  input: tuple val(outdir), path(primary_assembly), path(haplo_fasta), path(pacbio_reads)
  output: tuple path("primary_purged.fa"), path("haps_purged.fa"), path("*.log")
  script:
  """
  #! /usr/bin/env bash
  # === Inputs
  # genome_fasta=genome assembly file (cns_p_ctg.fasta)
  # haplo_fasta=haplotype genome (cns_h_ctg.fasta)
  # pacbio_reads=*subreads.fasta
  # === Outputs

  PROC=\$((`nproc`))

  ${samtools_app} faidx ${primary_assembly}
  
  for x in ${pacbio_reads.simpleName} ; do
    ${minimap2_app} \
      -xmap-pb \
      ${params.minimap2_params} \
      $primary_assembly \${x}.fasta \
      | $gzip_app -c - \
      > \${x}_p_mapping.paf.gz
  done
  
  ${pbcstat_app} *_p_mapping.paf.gz
  
  ${calcuts_app} PB.stat \
    > p_cutoffs \
    2> p_calcuts.log 
  
  ${split_fa_app} ${primary_assembly} > ${primary_assembly}.split
  
  ${minimap2_app} \
    -xasm5 \
    ${params.minimap2_params} \
    -DP \
    -t \${PROC} \
    ${primary_assembly}.split ${primary_assembly}.split \
    | gzip -c - \
    > ${primary_assembly}.split.self.paf.gz
  
  ${purge_dups_app} \
    -2 \
    -T p_cutoffs \
    -c PB.base.cov \
    ${primary_assembly}.split.self.paf.gz \
    > p_dups.bed \
    2> p_purge_dups.log
  
  ## -e flag to only remove haplotypic regions at the ends of contigs
  ${get_seqs_app} \
    -e \
    p_dups.bed ${primary_assembly} \
    -p primary

  
  echo "Purged primary, onto the alternate"
  
  ## do the same with the "alternate" haps after cating
  cat primary.hap.fa ${haplo_fasta} >  h_${haplo_fasta}
  
  ${samtools_app} faidx h_${haplo_fasta}
  
  for x in ${pacbio_reads.simpleName} ; do
    ${minimap2_app} \
      -xmap-pb \
      ${params.minimap2_params} \
      h_${haplo_fasta} \${x}.fasta \
      | $gzip_app -c - \
      > \${x}_h_mapping.paf.gz
  done
  
  ${pbcstat_app} *_h_mapping.paf.gz
  
  ${calcuts_app} PB.stat \
    > h_cutoffs \
    2> h_calcuts.log
  
  ${split_fa_app} h_${haplo_fasta} > h_${haplo_fasta}.split
  
  ${minimap2_app} \
    -xasm5 \
    ${params.minimap2_params} \
    -DP \
    -t \${PROC} \
    h_${haplo_fasta}.split h_${haplo_fasta}.split \
    | ${gzip_app} -c - \
    > h_${haplo_fasta}.split.self.paf.gz
  
  ${purge_dups_app} \
    -2 \
    -T h_cutoffs \
    -c PB.base.cov \
    h_${haplo_fasta}.split.self.paf.gz \
    > h_dups.bed \
    2> h_purge_dups.log

  ${get_seqs_app} \
    -e \
    h_dups.bed h_${haplo_fasta} \
    -p haps
  
  # rename to play nicely with nextflow simplename 
  mv primary.purged.fa primary_purged.fa
  mv haps.purged.fa haps_purged.fa
  """
  
  stub:
  """
  touch primary_purged.fa haps_purged.fa primary_hap.fa haps_hap.fa
  touch primary_purged.stats haps_purged.stats primary_hap.stats haps_hap.stats
  touch a.png a.log
  """
}

process PURGE_DUPS_CASE1 {
  publishDir "${params.outdir}/$outdir", mode:'copy'
  input: tuple val(outdir), path(primary_assembly), path(pacbio_reads)
  output: tuple path("primary_purged.fa"), path("haps_purged.fa"), path("*.log") //
  script:
  """
  #! /usr/bin/env bash
  # === Inputs
  # genome_fasta=genome assembly file (cns_p_ctg.fasta)
  # haplo_fasta=haplotype genome (cns_h_ctg.fasta)
  # pacbio_reads=*subreads.fasta
  # === Outputs
  ${samtools_app} faidx ${primary_assembly}
  
  for x in ${pacbio_reads.simpleName} ; do
    ${minimap2_app} \
    -xmap-pb \
    $primary_assembly \${x}.fasta \
    | $gzip_app -c - \
    > \${x}_p_mapping.paf.gz
  done
  
  ${pbcstat_app} *_p_mapping.paf.gz
  
  ${calcuts_app} PB.stat \
    > p_cutoffs \
    2> p_calcuts.log 
  
  ${split_fa_app} ${primary_assembly} > ${primary_assembly}.split
  
  ${minimap2_app} \
    -xasm5 \
    ${params.minimap2_params} \
    -DP \
    -t ${task.cpus} \
    ${primary_assembly}.split ${primary_assembly}.split \
    | gzip -c - \
    > ${primary_assembly}.split.self.paf.gz
  
  ${purge_dups_app} \
    -2 \
    -T p_cufoffs \
    -c PB.base.cov \
    ${primary_assembly}.split.self.paf.gz \
    > p_dups.bed \
    2> p_purge_dups.log
  
  ## -e flag to only remove haplotypic regions at the ends of contigs
  ${get_seqs_app} \
  -e \
  p_dups.bed ${primary_assembly} \
  -p primary
  
  echo "Purged primary, onto the alternate"
  
  ## do the same with the "alternate" haps after cating
  ${samtools_app} faidx primary.hap.fa
  
  for x in ${pacbio_reads.simpleName} ; do
    ${minimap2_app} \
      -xmap-pb \
      primary.hap.fa \${x}.fasta \
      | $gzip_app -c - \
      > \${x}_h_mapping.paf.gz
  done
  
  ${pbcstat_app} *_h_mapping.paf.gz
  
  ${calcuts_app} PB.stat \
    > h_cutoffs \
    2> h_calcuts.log
  
  ${split_fa_app} primary.hap.fa > primary.hap.split
  
  ${minimap2_app} \
    -xasm5 \
    -DP \
    -t ${task.cpus} \
    primary.hap.split primary.hap.split \
    | ${gzip_app} -c - \
    > primary.hap.split.self.paf.gz
  
  ${purge_dups_app} \
    -2 \
    -T h_cufoffs \
    -c PB.base.cov \
    primary.hap.split.self.paf.gz \
    > h_dups.bed \
    2> h_purge_dups.log
  
  if [[ -s h_dups.bed ]]; then
    ${get_seqs_app} \
      -e \
      h_dups.bed primary.hap.fa \
      -p haps
  else
    cat primary.hap.fa | sed 's/ REPEAT//g' > haps.purged.fa
  fi 
  
  ## rename to play nicely with nextflow simplename 
  mv primary.purged.fa primary_purged.fa
  mv haps.purged.fa haps_purged.fa
  """

  stub:
  """
  touch primary_purged.fa haps_purged.fa primary_hap.fa haps_hap.fa
  touch primary_purged.stats haps_purged.stats primary_hap.stats haps_hap.stats
  touch a.png a.log
  """
}

process PURGE_DUPS_TRIO {
  publishDir "${params.outdir}/$outdir", mode:'copy'
  input: tuple val(outdir), path(primary_assembly), path(pacbio_reads)
  output: tuple path("${primary_assembly.simpleName}_primary_purged.fa"), path("${primary_assembly.simpleName}_primary_hap.fa"), path("*.stats"), path("*.png"), path("*.log")
  script:
  template 'purge_dups_trios.sh'

  stub:
  """
  touch ${primary_assembly.simpleName}_primary_purged.fa ${primary_assembly.simpleName}_primary_hap.fa
  """
}


