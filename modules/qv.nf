#! /usr/bin/env nextflow

nextflow.enable.dsl=2

// create meryl database for merqury qv and merfin
process meryl_count {
  publishDir "${params.outdir}/00_Preprocess", mode: 'symlink'
  input: tuple val(k), path(illumina_read)
  output: path("*.meryl")
  script:
  """
  #! /usr/bin/env bash
  ${meryl_app} count k=${k} output ${illumina_read.simpleName}.meryl ${illumina_read}
  """

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
  """
  #! /usr/bin/env bash
  ${meryl_app} union-sum output illumina.meryl ${illumina_meryls}
  """

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
  """
  #! /usr/bin/env bash
  
  # Calculate histogram
  ${meryl_app} histogram ${illumina_meryl} > illumina.hist
  awk '\$1>5 {print}' illumina.hist | sort -k 2n | tail -n 1 | awk '{print \$1}' > peak.txt
  """

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
 \${merqury_sh} $illumina_db $assembly_fasta ${assembly_fasta.simpleName}

  # == Get single QV value and completeness stats
  cat ${assembly_fasta.simpleName}.qv | awk -F'\t' '{print \$4}' > merqury.qv
  cat ${assembly_fasta.simpleName}.completeness.stats | awk -F'\t' '{print \$5}' > comepleteness.qv
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
  """
  #! /usr/bin/env bash
  # Desc: Prints out length distibutions, GC, etc of each assembly, could be added to pipeline at the end
  # module load bbtools

  echo "Assmbly stats of $assembly_fasta  according to bbtools stats.sh"

  stats.sh in=$assembly_fasta out=${assembly_fasta.simpleName}.stats
  """

  stub:
  """
  touch bbstat_output.txt
  """
}
