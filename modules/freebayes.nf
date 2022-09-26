#! /usr/bin/env nextflow

nextflow.enable.dsl=2

include { create_windows; combineVCF; meryl_genome; merfin_polish; vcf_to_fasta } from './helper_functions.nf'

process align_shortreads {
  publishDir "${params.outdir}/${outdir}/bam", mode: 'symlink'
  input: tuple val(outdir), path(assembly_fasta), path(illumina_one), path(illumina_two)
  output: tuple val("$outdir"), path("*.bam"), path("*.bai")
  script:
  """
  #! /usr/bin/env bash
  PROC=\$(((`nproc`-1)*3/4+1))
  PROC2=\$(((`nproc`-1)*1/4+1))
  mkdir tmp
  ${bwamem2_app} index ${assembly_fasta}
  ${bwamem2_app} mem ${bwamem2_params} -t \$PROC ${assembly_fasta} ${illumina_one} ${illumina_two} | \
    ${samtools_app} sort -T tmp -m 8G --threads \$PROC2 - > ${illumina_one.simpleName}_aln.bam
  ${samtools_app} index -@ \${PROC} ${illumina_one.simpleName}_aln.bam
  """

  stub:
  """
  touch ${illumina_one.simpleName}_aln.bam
  touch ${illumina_one.simpleName}_aln.bam.bai
  """
}

process freebayes {
  errorStrategy { task.attempt < 4 ? 'retry' : 'terminate' }

  publishDir "${params.outdir}/${outdir}/vcf", mode: 'symlink'
  input: tuple val(outdir), path(illumina_bam), path(illumina_bai), path(assembly_fasta), path(assembly_fai), val(window)
  output: tuple val("$outdir"), path("*.vcf")
  script:
  """
  #! /usr/bin/env bash
  ${freebayes_app} \
    --region "${window}" \
    ${freebayes_params} \
    --bam ${illumina_bam} \
    --vcf ${illumina_bam.simpleName}_${window.replace(':','_').replace('|','_')}.vcf \
    --fasta-reference ${assembly_fasta}
  """

  stub:
  """
  touch ${illumina_bam.simpleName}_${window.replace(':','_').replace('|','_')}.vcf
  """
}

workflow FREEBAYES {
  take:
    outdir_ch
    asm_ch
    ill_ch
    peak_ch
    merylDB_ch
  main:
    win_ch = outdir_ch | combine(asm_ch) | create_windows | 
      map { n -> n.get(1) } | splitText() {it.trim() }
    fai_ch = create_windows.out | map { n -> n.get(0) }
    asm_meryl = outdir_ch | combine(channel.from(params.k)) | collect | combine(asm_ch) | meryl_genome

    new_asm_ch = outdir_ch | combine(asm_ch) | combine(ill_ch.collect()) | align_shortreads |
      combine(asm_ch) | combine(fai_ch) | combine(win_ch) |
      freebayes | groupTuple | combineVCF | combine(asm_ch) | 
      combine(asm_meryl) | combine(peak_ch) | combine(merylDB_ch) |  merfin_polish | combine(asm_ch) |
      vcf_to_fasta
  emit:
    new_asm_ch
}
