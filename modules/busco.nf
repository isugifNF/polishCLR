#! /usr/bin/env nextflow

nextflow.enable.dsl=2

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
  publishDir "${params.outdir}/${outdir}", mode: 'copy'
  scratch false

  input: tuple val(outdir), path(genomeFile)
  output: path("${genomeFile.simpleName}/*")

  script:
  """
  #! /usr/bin/env bash
  PROC=\$((`nproc`))

  if [ -s ${genomeFile} ]; then

  cat ${genomeFile} | tr '|' '_' > ${genomeFile.simpleName}_fixheaders.fna

  ${busco_app} \
    -o ${genomeFile.simpleName} \
    -i ${genomeFile.simpleName}_fixheaders.fna \
    ${busco_params} \
    -c \${PROC}

else
  mkdir ${genomeFile.simpleName}
  echo '${genomeFile.simpleName} is empty, this is a stub' > ${genomeFile.simpleName}/empty.txt

fi
  """

  stub:
  """
  mkdir ${genomeFile.simpleName}
  touch ${genomeFile.simpleName}/hey.txt
  """
}
