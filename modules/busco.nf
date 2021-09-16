#! /usr/bin/env nextflow

nextflow.enable.dsl=2

MERQURY='' // define to fix clash with other conda

// process BUSCO_setup {
//   publishDir "${params.outdir}/03c_BUSCO", mode: 'copy'
//   output:tuple path("config"), path("Busco_version.txt")
// 
//   script:
//   """
//   #! /usr/bin/env bash
//   cp -r /augustus/config .
//   ${busco_app} --version > Busco_version.txt
//   """
// }

process BUSCO {
  publishDir "${params.outdir}/03c_BUSCO", mode: 'copy'
  scratch false

  input: path(genomeFile)
  output: path("${genomeFile.simpleName}/*")
  
  script:
  """
  #! /usr/bin/env bash
  PROC=\$((`nproc`))
  cat ${genomeFile} | tr '|' '_' > ${genomeFile.simpleName}_fixheaders.fna
  ${busco_app} \
    -o ${genomeFile.simpleName} \
    -i ${genomeFile.simpleName}_fixheaders.fna \
    -l ${busco_lineage} \
    -m genome \
    -c \${PROC} \
    -f
  """
}