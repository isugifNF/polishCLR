#! /usr/bin/env nextflow

nextflow.enable.dsl=2

// create meryl database for merqury qv and merfin
process meryl_count {
  publishDir "${params.outdir}/00_Preprocess", mode: 'symlink'
  input: tuple val(k), path(illumina_read)
  output: path("*.meryl")
  script:
  template 'meryl_count.sh'

  stub:
  """
  touch ${illumina_read.simpleName}.meryl
  """
}

process meryl_union {
  publishDir "${params.outdir}/00_Preprocess", mode: 'copy'
  input: path(illumina_meryls)
  output: path("illumina.meryl")
  script:
  template 'meryl_union.sh'

  stub:
  """
  touch illumina.meryl
  """
}

// TODO: combine the above in a mk_meryldb workflow

// calculate illumina peak for merfin
process meryl_peak {
  publishDir "${params.outdir}/00_Preprocess", mode: 'symlink'
  input: path(illumina_meryl)
  output: tuple path("peak.txt"), path("illumina.hist")
  script:
  template 'meryl_peak.sh'

  stub:
  """
  touch peak.txt illumina.hist
  echo "72" >> peak.txt
  """
}

// 01 Merqury QV value
process MerquryQV {
  publishDir "${params.outdir}/${outdir}/${assembly_fasta.simpleName}/MerquryQV", mode: 'copy'
  publishDir "${params.outdir}/${outdir}/${assembly_fasta.simpleName}", mode: 'copy', pattern: "merqury.qv"

  input: tuple val(outdir), path(illumina_db), path(assembly_fasta)
  output: path("*")
  script:
  template 'merquryqv.sh'

  stub:
  """
  touch merqury.qv
  """
}

// 01 bbstat: Length distribtions of initial assembly
process bbstat {
  publishDir "${params.outdir}/${outdir}/${assembly_fasta.simpleName}/bbstat", mode: 'copy'
  input: tuple val(outdir), path(assembly_fasta)
  output: path("*")
  script:
  template 'bbstats.sh'

  stub:
  """
  touch bbstat_output.txt
  """
}