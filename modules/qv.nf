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
  publishDir "${params.outdir}/${outdir}/MerquryQV", mode: 'copy'
  publishDir "${params.outdir}/${outdir}/", mode: 'copy', pattern: "merqury.qv"

  input: tuple val(outdir), path(illumina_db), path(assembly_fasta)
  output: path("*")
  script:
  """
  #! /usr/bin/env bash

  merqury_sh="$params.merqury_sh"

  printf "======================================================= \n"
  printf "uniq kmers in asm | kmers in both asm and reads | QV | Error rate \n"
  printf "======================================================= \n"
  \${merqury_sh} $illumina_db $assembly_fasta ${assembly_fasta.simpleName}
  printf "======================================================= \n" > ${assembly_fasta.simpleName}_qv.txt
  printf "uniq kmers in asm | kmers in both asm and reads | QV | Error rate \n" >> ${assembly_fasta.simpleName}_qv.txt
  printf "======================================================= \n" >> ${assembly_fasta.simpleName}_qv.txt
  cat  ${assembly_fasta.simpleName}.qv >> ${assembly_fasta.simpleName}_qv.txt

  # == Get single QV value
  cat ${assembly_fasta.simpleName}.qv | awk -F'\t' '{print \$4}' > merqury.qv
  """

  stub:
  """
  touch merqury.qv other.png
  """
}

// 01 bbstat: Length distribtions of initial assembly
process bbstat {
  publishDir "${params.outdir}/${outdir}", mode: 'copy'
  input: tuple val(outdir), path(assembly_fasta)
  output: path("*")
  script:
  template 'bbstats.sh'

  stub:
  """
  touch bbstat_output.txt
  """
}