#! /usr/bin/env nextflow

nextflow.enable.dsl=2

include { create_windows;
          combineVCF;
          meryl_genome;
          merfin_polish;
          vcf_to_fasta } from './helper_functions.nf'

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
    win_ch = outdir_ch
      | combine(asm_ch)
      | create_windows
      | map { n -> n.get(1) }
      | splitText() {it.trim() }

    fai_ch = create_windows.out
      | map { n -> n.get(0) }

    newasm_ch = outdir_ch
      | combine(asm_ch)
      | pbmm2_index
      | combine(pac_ch)
      | pbmm2_align
      | combine(asm_ch)
      | combine(fai_ch)
      | combine(win_ch)
      | gcpp_arrow
      | map { n -> [ n.get(0), n.get(1) ] }
      | groupTuple
      | merge_consensus

  emit:
    newasm_ch
}

// 2nd Arrow run with merfin

process reshape_arrow {
  publishDir "${params.outdir}/${outdir}/merfin", mode: 'symlink'
  input: tuple val(outdir), path(vcf)
  output: tuple val("${outdir}"), path("*.reshaped.vcf.gz")
  script:
  """
  #! /usr/bin/env bash

  # create the extra_header.vcf
  cat << EOF > extra_header.vcf 
  ##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
  ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
  ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
  EOF
  
  #  reshape_arrow.sh
  # https://github.com/arangrhie/merfin/blob/master/scripts/reformat_arrow/reshape_arrow.sh
  output=${vcf.baseName} # assumes input.vcf
  echo \${output}
  grep -v "#" \${output}.vcf | sed 's/,/;/g' > \${output}.temp.reshaped.vcf
  ${bcftools_app} view -h \${output}.vcf > \${output}.temp.reshaped.header.vcf
  cat \${output}.temp.reshaped.header.vcf \${output}.temp.reshaped.vcf > \${output}.temp.reshaped.combined.vcf
  rm \${output}.temp.reshaped.header.vcf \${output}.temp.reshaped.vcf
  ${bcftools_app} annotate -h extra_header.vcf \${output}.temp.reshaped.combined.vcf > \${output}.temp.reshaped.vcf
  ${bcftools_app} view -h \${output}.temp.reshaped.vcf | sed 's/\tINFO/\tINFO\tFORMAT\tIND/g' > \${output}.reshaped.vcf
  rm \${output}.temp.reshaped.vcf 
  ${bcftools_app} view -H \${output}.temp.reshaped.combined.vcf | awk -F"\t" -v OFS="\t" '{gsub(/DP=/,".\tGT:DP\t1/1:",\$8);print \$0}' >> \${output}.reshaped.vcf

  # https://github.com/arangrhie/merfin/wiki/Best-practices-for-Merfin#preparing-input-vcf
  ${bcftools_app} view -Oz -i '(GT="AA" || GT="Aa")' \${output}.reshaped.vcf > out.vcf.gz
  ${bcftools_app} view out.vcf.gz -Oz > \${output}.reshaped.vcf.gz
  rm \${output}.reshaped.vcf 
  rm \${output}.temp.reshaped.combined.vcf
  """

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
    win_ch = outdir_ch
      | combine(asm_ch)
      | create_windows
      | map { n -> n.get(1) }
      | splitText() {it.trim() }

    fai_ch = create_windows.out
      | map { n -> n.get(0) }

    arrow_run_ch = outdir_ch
      | combine(asm_ch)
      | pbmm2_index
      | combine(pac_ch)
      | pbmm2_align
      | combine(asm_ch)
      | combine(fai_ch)
      | combine(win_ch)
      | gcpp_arrow

    if ( params.same_specimen ) {
      /* create a genome meryl db */
      asm_meryl = outdir_ch
        | combine(channel.of(params.k))
        | combine(asm_ch)
        | meryl_genome

      /* prepare and run merfin polish */
      newasm_ch = arrow_run_ch
        | map { n -> [ n.get(0), n.get(2) ] }
        | groupTuple
        | combineVCF
        | reshape_arrow
        | combine(asm_ch)
        | combine(asm_meryl)
        | combine(peak_ch)
        | combine(merylDB_ch)
        | merfin_polish
        | combine(asm_ch)
        | vcf_to_fasta
    } else {
      newasm_ch = arrow_run_ch
        | map { n -> [ n.get(0), n.get(1) ] }
        | groupTuple
        | merge_consensus
    }

  emit:
    newasm_ch
}
