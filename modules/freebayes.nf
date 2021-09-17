#! /usr/bin/env nextflow

nextflow.enable.dsl=2

process align_shortreads {
  publishDir "${params.outdir}/${outdir}/bam", mode: 'symlink'
  input: tuple val(outdir), path(assembly_fasta), path(illumina_one), path(illumina_two)
  output: tuple val("$outdir"), path("*.bam"), path("*.bai")
  script:
  template 'align_shortreads.sh'
}

process freebayes {
  publishDir "${params.outdir}/${outdir}/vcf", mode: 'symlink'
  input: tuple val(outdir), path(illumina_bam), path(illumina_bai), path(assembly_fasta), path(assembly_fai), val(window)
  output: tuple val("$outdir"), path("*.vcf")
  script:
  template 'freebayes.sh'
}

process combineVCF {
  publishDir "${params.outdir}/${outdir}", mode: 'symlink'
  input: tuple val(outdir), path(vcfs)
  output: tuple val("${outdir}"), path("${i}_consensus.vcf")
  script:
  template 'combineVCF.sh'
}

process meryl_genome_fb {
  publishDir "${params.outdir}/${outdir}/merfin", mode: 'symlink'
  input: tuple val(outdir), val(k), path(illumina_read)
  output: path("*.meryl")
  script:
  template 'meryl_count.sh'
}

process merfin_polish {
  publishDir "${params.outdir}/${outdir}/merfin", mode: 'symlink'
  input: tuple val(outdir), path(vcf), path(genome_fasta), path(genome_meryl), val(peak), path(meryldb)
  output: tuple val("$outdir"), path("*merfin.polish.vcf")
  script:
  template 'merfin_polish.sh'
}

process vcf_to_fasta {
  publishDir "${params.outdir}/${outdir}", mode: 'symlink'
  input: tuple val(outdir), path(vcf), path(genome_fasta)
  output: path("${outdir}_consensus.fasta")
  script:
  template 'vcf_to_fasta.sh'
}

workflow FREEBAYES {
  take:
    outdir_ch // 06 or 08_FreeBayesPolish
    asm_ch
    ill_ch
    peak_ch
    merylDB_ch
  main:
    win_ch = outdir_ch | combine(asm_ch) | create_windows | 
      map { n -> n.get(1) } | splitText() {it.trim() }
    fai_ch = create_windows.out | map { n -> n.get(0) }
    asm_meryl = outdir_ch | combine(channel.from(params.k)) | collect | combine(asm_ch) | meryl_genome_fb

    new_asm_ch = outdir_ch | combine(asm_ch) | combine(ill_ch.collect()) | align_shortreads |
      combine(asm_ch) | combine(fai_ch) | combine(win_ch) |
      freebayes | groupTuple | combineVCF | combine(asm_ch) | 
      combine(asm_meryl) | combine(peak_ch) | combine(merylDB_ch) |  merfin_polish | combine(asm_ch) |
      vcf_to_fasta
  emit:
    new_asm_ch
}