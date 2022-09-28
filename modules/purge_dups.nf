#! /usr/bin/env nextflow

nextflow.enable.dsl=2

process PURGE_DUPS {
  publishDir "${params.outdir}/$outdir", mode:'copy'
  input: tuple val(outdir), path(primary_assembly), path(haplo_fasta), path(pacbio_reads)
  output: tuple path("primary_purged.fa"), path("haps_purged.fa"), path("*.log") //
  script:
  template 'purge_dups.sh'

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
  template 'purge_dups_case1.sh'

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


