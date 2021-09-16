#! /usr/bin/env nextflow

nextflow.enable.dsl=2

// create meryl database for merqury qv and merfin
process meryl_count {
  publishDir "${params.outdir}/00_Preprocess", mode: 'symlink'
  input: tuple val(k), path(illumina_read)
  output: path("*.meryl")
  script:
  template 'meryl_count.sh'
}

process meryl_union {
  publishDir "${params.outdir}/00_Preprocess", mode: 'copy'
  input: path(illumina_meryls)
  output: path("illumina.meryl")
  script:
  template 'meryl_union.sh'
}

// calculate illumina peak for merfin
process meryl_peak {
  publishDir "${params.outdir}/00_Preprocess", mode: 'symlink'
  input: path(illumina_meryl)
  output: tuple path("peak.txt"), path("illumina.hist")
  script:
  template 'meryl_peak.sh'
}

// 01 Merqury QV value
process MerquryQV {
  publishDir "${params.outdir}/01_QV/MerquryQV", mode: 'copy'
  publishDir "${params.outdir}/01_QV", mode: 'copy', pattern: "merqury.qv"
  input: tuple path(illumina_db), path(assembly_fasta)
  output: path("*")
  script:
  template 'merquryqv.sh'
}

// 01 bbstat: Length distribtions of initial assembly
process bbstat {
  publishDir "${params.outdir}/01_QV/bbstat", mode: 'copy'
  input:  path(assembly_fasta)
  output: path("*")
  script:
  template 'bbstats.sh'
}