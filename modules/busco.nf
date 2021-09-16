#! /usr/bin/env nextflow

nextflow.enable.dsl=2

// Temporary solution till mapping to path
busco_container = 'ezlabgva/busco:v5.1.2_cv1'

// this setup is required because BUSCO runs Augustus that requires writing to the config/species folder.  
// So this folder must be bound outside of the container and therefore needs to be copied outside the container first.
// Leaving in container code temporarily
process BUSCO_setup {
  publishDir "${params.outdir}/03c_BUSCO", mode: 'copy'
  container = "$busco_container"

  output:tuple path("config"), path("Busco_version.txt")

  script:
  """
  #! /usr/bin/env bash
  cp -r /augustus/config .
  echo "$busco_container" > Busco_version.txt
  """
}

process BUSCO {
  publishDir "${params.outdir}/03c_BUSCO", mode: 'copy'
  container = "$busco_container"
  scratch false

  input: tuple path(genomeFile), path(config)
  output: tuple  path("${genomeFile.simpleName}/*")
  
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