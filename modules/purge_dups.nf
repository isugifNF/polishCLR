#! /usr/bin/env nextflow

nextflow.enable.dsl=2

process PURGE_DUPS {
  publishDir "${params.outdir}/03b_PurgeDups", mode:'copy'
  input: tuple path(primary_assembly), path(haplo_fasta), path(pacbio_reads)
  output: tuple path("primary_purged.fa"), path("haps_purged.fa"), path("primary_hap.fa"), path("haps_hap.fa") //, path("h_${haplo_fasta}") //, path("*.stats")
  script:
  template 'purge_dups.sh'
}