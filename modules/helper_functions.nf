#! /usr/bin/env nextflow

nextflow.enable.dsl=2

// 00 Preprocess: Unzip any bz2 files
process bz_to_gz {
  publishDir "${params.outdir}/00_Preprocess", mode: 'copy'
  input:tuple val(readname), path(illumina_reads)
  output: tuple val(readname), path("*.gz")
  script:
  template 'bz_to_gz.sh'
}

// concat genome and mito together
process MERGE_FILE {
  publishDir "${params.outdir}/00_Preprocess", mode: 'copy'
  input: tuple path(primary_assembly), path(alternate_assembly), path(mito_assembly)
  output: path("${primary_assembly.simpleName}_merged.fasta")
  script:
  template 'merge_file.sh'
}

process MERGE_FILE_TRIO {
  publishDir "${params.outdir}/00_Preprocess", mode: 'copy'
  input: tuple path(primary_assembly), path(mito_assembly)
  output: path("${primary_assembly.simpleName}_merged.fasta")
  script:
  template 'merge_file_trio.sh'
}

process SPLIT_FILE {
  input: path(genome_fasta)
  output: tuple path("p_${genome_fasta}"), path("a_${genome_fasta}"), path("m_${genome_fasta}")
  script:
  template 'split_file.sh'
}