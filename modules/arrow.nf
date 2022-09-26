#! /usr/bin/env nextflow

nextflow.enable.dsl=2

include { create_windows; combineVCF; meryl_genome; merfin_polish; vcf_to_fasta } from './helper_functions.nf'

process pbmm2_index {
  publishDir "${params.outdir}/${outdir}", mode: 'symlink'
  input: tuple val(outdir), path(assembly_fasta)
  output: tuple val("$outdir"), path("$assembly_fasta"), path("*.mmi")
  script:
  """
  #! /usr/bin/env bash
  ${pbmm2_app} index ${assembly_fasta} ${assembly_fasta}.mmi
  """

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
  """
  #! /usr/bin/env bash
  PROC=\$(((`nproc`-1)*3/4+1))
  PROC2=\$(((`nproc`-1)*1/4+1))
  mkdir tmp

  # for multiple pacbio subread files
  ls ${pacbio_read} > bam.fofn

  ${pbmm2_app} align ${pbmm2_params} -j \$PROC ${assembly_fasta} bam.fofn | \
    ${samtools_app} sort -T tmp -m 8G --threads \$PROC2 - > ${assembly_fasta.simpleName}_aln.bam
  ${samtools_app} index -@ \${PROC} ${assembly_fasta.simpleName}_aln.bam
  """

  stub:
  """
  touch ${assembly_fasta.simpleName}_${pacbio_read.simpleName}_aln.bam
  touch ${assembly_fasta.simpleName}_${pacbio_read.simpleName}_aln.bam.bai
  """
}

process gcpp_arrow {
  errorStrategy { task.attempt < 4 ? 'retry' : 'terminate' }

  publishDir "${params.outdir}/${outdir}/gccpruns", mode: 'symlink'
  input: tuple val(outdir), path(pacbio_bam), path(pacbio_bai),  path(assembly_fasta), path(assembly_fai), val(window)
  output: tuple val("$outdir"), path("*.fasta"), path("*.vcf")
  script:
  """
  #! /usr/bin/env bash
  PROC=\$(((`nproc`-1)*3/4+1))
  ${gcpp_app} --algorithm=arrow \
    ${gcpp_params} \
    -j \${PROC} -w "$window" \
    -r ${assembly_fasta} ${pacbio_bam} \
    -o ${assembly_fasta.simpleName}_${window.replace(':','_').replace('|','_')}.vcf,${assembly_fasta.simpleName}_${window.replace(':','_').replace('|','_')}.fasta
  """

  stub:
  """
  touch ${assembly_fasta.simpleName}_${window.replace(':','_').replace('|','_')}.vcf
  touch ${assembly_fasta.simpleName}_${window.replace(':','_').replace('|','_')}.fasta
  """
}

process merge_consensus {
  publishDir "${params.outdir}/${outdir}", mode: 'copy'
  input: tuple val(outdir), path(windows_fasta)
  output: path("*_consensus.fasta")
  script:
  """
  #! /usr/bin/env bash
  OUTNAME=`echo "$outdir" | sed 's:/:_:g'`
  cat ${windows_fasta} > \${OUTNAME}_consensus.fasta
  """

  stub:
  """
  OUTNAME=`echo "$outdir" | sed 's:/:_:g'`
  touch \${OUTNAME}_consensus.fasta
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
      combine(asm_ch) | combine(fai_ch) | combine(win_ch) | gcpp_arrow | 
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
      combine(asm_ch) | combine(fai_ch) | combine(win_ch) | gcpp_arrow 
    if (params.same_specimen) {
      /* create a genome meryl db */
      asm_meryl = outdir_ch | combine(channel.of(params.k)) | combine(asm_ch) | meryl_genome 

      /* prepare and run merfin polish */
      newasm_ch = arrow_run_ch | map { n -> [ n.get(0), n.get(2) ] } | groupTuple |
        combineVCF | reshape_arrow | combine(asm_ch) |  combine(asm_meryl) |combine(peak_ch) |
        combine(merylDB_ch) | merfin_polish | combine(asm_ch) | vcf_to_fasta
    } else {
      newasm_ch = arrow_run_ch | map { n -> [ n.get(0), n.get(1) ] } | groupTuple |
        merge_consensus
    }
  
  emit:
    newasm_ch
}
