#! /usr/bin/env nextflow

nextflow.enable.dsl=2

process create_windows {
  publishDir "${params.outdir}/0${i}_ArrowPolish", mode: 'symlink'
  input: tuple val(i), path(assembly_fasta)
  output: tuple path("*.fai"), path("win.txt")
  script:
  template 'create_windows.sh'
}

process pbmm2_index {
  publishDir "${params.outdir}/0${i}_ArrowPolish", mode: 'symlink'
  input: tuple val(i), path(assembly_fasta)
  output: tuple val("$i"), path("$assembly_fasta"), path("*.mmi")
  script:
  template 'pbmm2_index.sh'
}

process pbmm2_align {
  publishDir "${params.outdir}/0${i}_ArrowPolish", mode: 'symlink'
  input:tuple val(i), path(assembly_fasta), path(assembly_mmi), path(pacbio_read)
  output: tuple val("$i"), path("*.bam"), path("*.bai")
  script:
  template 'pbmm2_align.sh'
}

process gcc_Arrow {
  publishDir "${params.outdir}/0${i}_ArrowPolish/gccruns", mode: 'symlink'
  input: tuple val(i), path(pacbio_bam), path(pacbio_bai),  path(assembly_fasta), path(assembly_fai), val(window)
  output: tuple val("$i"), path("*.fasta"), path("*.vcf")
  script:
  template 'gcc_arrow.sh'
}

process merge_consensus {
  publishDir "${params.outdir}/0${i}_ArrowPolish", mode: 'copy'
  input: tuple val(i), path(windows_fasta)
  output: path("${i}_consensus.fasta")
  script:
  template 'merge_consensus.sh'
}

workflow ARROW {
  take:
    asm_ch
    pac_ch
  main:
    win_ch = channel.of("2") | combine(asm_ch) | create_windows | 
      map { n -> n.get(1) } | splitText() {it.trim() }
    fai_ch = create_windows.out | map { n -> n.get(0) }

    newasm_ch = channel.of("2") | combine(asm_ch) | pbmm2_index | combine(pac_ch) | pbmm2_align |
      combine(asm_ch) | combine(fai_ch) | combine(win_ch) | gcc_Arrow | 
      map { n -> [ n.get(0), n.get(1) ] } | groupTuple | 
      merge_consensus
  
  emit:
    newasm_ch
}

// 2nd Arrow run with merfin
process meryl_genome {
  publishDir "${params.outdir}/04_ArrowPolish/merfin", mode: 'symlink'
  input: tuple val(k), path(illumina_read)
  output: path("*.meryl")
  script:
  template 'meryl_count.sh'
}

process combineVCF_arrow {
  publishDir "${params.outdir}/0${i}_ArrowPolish", mode: 'symlink'
  input: tuple val(i), path(vcfs)
  output: tuple val(i), path("${i}_consensus.vcf")
  script:
  template 'combineVCF.sh'
}

process reshape_arrow {
  publishDir "${params.outdir}/0${i}_ArrowPolish/merfin", mode: 'symlink'
  input: tuple val(i), path(vcf)
  output: tuple val(i), path("*.reshaped.vcf.gz")
  script:
  template 'reshape_arrow.sh'
}

process merfin_polish_arrow {
  publishDir "${params.outdir}/0${i}_ArrowPolish/merfin", mode: 'symlink'
  input: tuple val(i), path(vcf), path(genome_fasta), path(genome_meryl), val(peak), path(meryldb)
  output: tuple val("$i"), path("*merfin.polish.vcf")
  script:
  template 'merfin_polish.sh'
}

process vcf_to_fasta_arrow {
  publishDir "${params.outdir}/0${i}_ArrowPolish", mode: 'symlink'
  input: tuple val(i), path(vcf), path(genome_fasta)
  output: path("${i}_consensus.fasta")
  script:
  template 'vcf_to_fasta.sh'
}

workflow ARROW_MERFIN {
  take:
    asm_ch
    pac_ch
    peak_ch
    merylDB_ch
    
  main:
    win_ch = channel.of("4") | combine(asm_ch) | create_windows | 
      map { n -> n.get(1) } | splitText() {it.trim() }
    fai_ch = create_windows.out | map { n -> n.get(0) }

    arrow_run_ch = channel.of("4") | combine(asm_ch) | pbmm2_index | combine(pac_ch) | pbmm2_align |
      combine(asm_ch) | combine(fai_ch) | combine(win_ch) | gcc_Arrow 

    if (params.same_specimen) {
      /* create a genome meryl db */
      asm_meryl = channel.of(params.k) | combine(asm_ch) | meryl_genome

      /* prepare and run merfin polish */
      newasm_ch = arrow_run_ch | map { n -> [ n.get(0), n.get(2) ] } | groupTuple |
        combineVCF_arrow | reshape_arrow | combine(asm_ch) | combine(asm_meryl) | combine(peak_ch) |
        combine(merylDB_ch) | merfin_polish_arrow | combine(asm_ch) | vcf_to_fasta_arrow
    } else {
      newasm_ch = arrow_run_ch | map { n -> [ n.get(0), n.get(1) ] } | groupTuple |
        merge_consensus
    }
  
  emit:
    newasm_ch
}