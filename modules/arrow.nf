#! /usr/bin/env nextflow

nextflow.enable.dsl=2

include { create_windows; combineVCF; meryl_genome; merfin_polish; vcf_to_fasta } from './helper_functions.nf'

process pbmm2_index {
  publishDir "${params.outdir}/${outdir}", mode: 'symlink'
  input: tuple val(outdir), path(assembly_fasta)
  output: tuple val("$outdir"), path("$assembly_fasta"), path("*.mmi")
  script:
  template 'pbmm2_index.sh'

  stub:
  """
  touch ${assembly_fasta}
  touch ${assembly_fasta}.mmi
  """
}

process pbmm2_align {
  publishDir "${params.outdir}/${outdir}", mode: 'symlink'
  input:tuple val(outdir), path(assembly_fasta), path(assembly_mmi), path(pacbio_read)
  output: tuple val("$outdir"), path("*.bam"), path("*.bai")
  script:
  template 'pbmm2_align.sh'

  stub:
  """
  touch ${pacbio_read.simpleName}_aln.bam
  touch ${pacbio_read.simpleName}_aln.bai
  """
}

process gcc_Arrow {
  publishDir "${params.outdir}/${outdir}/gccruns", mode: 'symlink'
  input: tuple val(outdir), path(pacbio_bam), path(pacbio_bai),  path(assembly_fasta), path(assembly_fai), val(window)
  output: tuple val("$outdir"), path("*.fasta"), path("*.vcf")
  script:
  template 'gcc_arrow.sh'

  stub:
  """
  touch ${window.replace(':','_').replace('|','_')}.vcf
  touch ${window.replace(':','_').replace('|','_')}.fasta
  """
}

process merge_consensus {
  publishDir "${params.outdir}/${outdir}", mode: 'copy'
  input: tuple val(outdir), path(windows_fasta)
  output: path("${outdir}_consensus.fasta")
  script:
  template 'merge_consensus.sh'

  stub:
  """
  touch ${outdir}_consensus.fasta
  """
}

workflow ARROW {
  take:
    outdir_ch // 02_ArrowPolish
    asm_ch
    pac_ch
  main:
    win_ch = outdir_ch | combine(asm_ch) | create_windows | 
      map { n -> n.get(1) } | splitText() {it.trim() }
    fai_ch = create_windows.out | map { n -> n.get(0) }

    newasm_ch = outdir_ch | combine(asm_ch) | pbmm2_index | combine(pac_ch) | pbmm2_align |
      combine(asm_ch) | combine(fai_ch) | combine(win_ch) | gcc_Arrow | 
      map { n -> [ n.get(0), n.get(1) ] } | groupTuple | 
      merge_consensus
  
  emit:
    newasm_ch
}

// 2nd Arrow run with merfin

process reshape_arrow {
  publishDir "${params.outdir}/${outdir}/merfin", mode: 'symlink'
  input: tuple val(outdir), path(vcf)
  output: tuple val("${outdir}"), path("*.reshaped.vcf.gz")
  script:
  template 'reshape_arrow.sh'

  stub:
  """
  touch ${vcf.baseName}.reshaped.vcf.gz
  """
}

workflow ARROW_MERFIN {
  take:
    outdir_ch   // 04_ArrowPolish
    asm_ch
    pac_ch
    peak_ch
    merylDB_ch
    
  main:
    win_ch = outdir_ch | combine(asm_ch) | create_windows | 
      map { n -> n.get(1) } | splitText() {it.trim() }
    fai_ch = create_windows.out | map { n -> n.get(0) }

    arrow_run_ch = outdir_ch | combine(asm_ch) | pbmm2_index | combine(pac_ch) | pbmm2_align |
      combine(asm_ch) | combine(fai_ch) | combine(win_ch) | gcc_Arrow 

    if (params.same_specimen) {
      /* create a genome meryl db */
      asm_meryl = outdir_ch | combine(channel.of(params.k)) | combine(asm_ch) | meryl_genome

      /* prepare and run merfin polish */
      newasm_ch = arrow_run_ch | map { n -> [ n.get(0), n.get(2) ] } | groupTuple |
        combineVCF | reshape_arrow | combine(asm_ch) | combine(asm_meryl) | combine(peak_ch) |
        combine(merylDB_ch) | merfin_polish | combine(asm_ch) | vcf_to_fasta
    } else {
      newasm_ch = arrow_run_ch | map { n -> [ n.get(0), n.get(1) ] } | groupTuple |
        merge_consensus
    }
  
  emit:
    newasm_ch
}