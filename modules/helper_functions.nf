#! /usr/bin/env nextflow

nextflow.enable.dsl=2

// 00 Preprocess: Unzip any bz2 files
process bz_to_gz {
  publishDir "${params.outdir}/00_Preprocess", mode: 'copy'
  input:tuple val(readname), path(illumina_reads)
  output: tuple val(readname), path("*.gz")
  script:
  """
  #! /usr/bin/env bash

  PROC=\$(((`nproc`-1)/2+1))
  ${parallel_app} ${parallel_params} "${bzcat_app} {1} | ${pigz_app} -p \${PROC} > {1/.}.gz" ::: *.bz2
  """

  stub:
  """
  touch ${illumina_reads.get(0)}.gz
  touch ${illumina_reads.get(1)}.gz
  """
}

// concat genome and mito together
process MERGE_FILE {
  publishDir "${params.outdir}/00_Preprocess", mode: 'copy'
  input: tuple path(primary_assembly), path(alternate_assembly), path(mito_assembly)
  output: path("${primary_assembly.simpleName}_merged.fasta")
  script:
  """
  #! /usr/bin/env bash

  # === Inputs
  # primary_assembly = p_ctg.fasta     # From FALCON or FALCON-Unzip
  # alternate_assembly = a_ctg.fasta
  # mito_assembly = mt.fasta           # From vgpMito pipeline

  # === Outputs
  # ${primary_assembly.simpleName}_merged.fasta

  cat ${primary_assembly} | sed 's/>/>pri_/g' > ${primary_assembly.simpleName}_temp.fasta
  echo "" >> ${primary_assembly.simpleName}_temp.fasta
  cat ${mito_assembly} | sed 's/>/>mit_/g' >> ${primary_assembly.simpleName}_temp.fasta
  echo "" >> ${primary_assembly.simpleName}_temp.fasta
  cat ${alternate_assembly} | sed 's/>/>alt_/g' >> ${primary_assembly.simpleName}_temp.fasta
  cat ${primary_assembly.simpleName}_temp.fasta | grep -v "^\$" > ${primary_assembly.simpleName}_merged.fasta
  """

  stub:
  """
  touch ${primary_assembly.simpleName}_merged.fasta
  """

}

// concat genome and mito together
process MERGE_FILE_CASE1 {
  publishDir "${params.outdir}/00_Preprocess", mode: 'copy'
  input: tuple path(primary_assembly), path(mito_assembly)
  output: path("${primary_assembly.simpleName}_merged.fasta")
  script:
  """
  #! /usr/bin/env bash

  # === Inputs
  # primary_assembly = p_ctg.fasta     # From CANU, similar to Falcon but without alternative assembly
  # mito_assembly = mt.fasta           # From vgpMito pipeline

  # === Outputs
  # ${primary_assembly.simpleName}_merged.fasta

  cat ${primary_assembly} | sed 's/>/>pri_/g' > ${primary_assembly.simpleName}_temp.fasta
  echo "" >> ${primary_assembly.simpleName}_temp.fasta
  cat ${mito_assembly} | sed 's/>/>mit_/g' >> ${primary_assembly.simpleName}_temp.fasta
  cat ${primary_assembly.simpleName}_temp.fasta | grep -v "^\$" > ${primary_assembly.simpleName}_merged.fasta
  """

  stub:
  """
  touch ${primary_assembly.simpleName}_merged.fasta
  """

}

process SPLIT_FILE_p {
  publishDir "${params.outdir}/${outdir}", mode:'copy'
  input: path(genome_fasta)
  output: tuple path("pat_${genome_fasta}"), path("mit_${genome_fasta}")
  script:
  """
  #! /usr/bin/env bash

  # === Inputs
  # genome_fasta = primary_assembly_merged.fasta
  # === Outputs
  # p_${genome_fasta}     # primary assembly
  # a_${genome_fasta}     # alternative assembly
  # m_${genome_fasta}     # mitochondrial assembly

  ${samtools_app} faidx ${genome_fasta}
  grep ">pri_" ${genome_fasta} | cut -f1 | sed 's/>//g' > pri.list
  ${samtools_app} faidx -r pri.list ${genome_fasta} > pat_${genome_fasta}

  grep ">mit_" ${genome_fasta} | cut -f1 | sed 's/>//g' > mit.list
  ${samtools_app} faidx -r mit.list ${genome_fasta} > mit_${genome_fasta}
  """

  stub:
  """
  touch pat_${genome_fasta} mit_${genome_fasta}
  """
}

process SPLIT_FILE_m {
  publishDir "${params.outdir}/${outdir}", mode:'copy'
  input: path(genome_fasta)
  output: tuple path("mat_${genome_fasta}"), path("mit_${genome_fasta}")
  script:
  """
  #! /usr/bin/env bash

  # === Inputs
  # genome_fasta = primary_assembly_merged.fasta
  # === Outputs
  # p_${genome_fasta}     # primary assembly
  # a_${genome_fasta}     # alternative assembly
  # m_${genome_fasta}     # mitochondrial assembly

  ${samtools_app} faidx ${genome_fasta}
  grep ">pri_" ${genome_fasta} | cut -f1 | sed 's/>//g' > pri.list
  ${samtools_app} faidx -r pri.list ${genome_fasta} > mat_${genome_fasta}

  grep ">mit_" ${genome_fasta} | cut -f1 | sed 's/>//g' > mit.list
  ${samtools_app} faidx -r mit.list ${genome_fasta} > mit_${genome_fasta}
  """

  stub:
  """
  touch mat_${genome_fasta} mit_${genome_fasta}
  """
}


process SPLIT_FILE_CASE1 {
  publishDir "${params.outdir}/${outdir}", mode:'copy'

  input:tuple val(outdir), path(genome_fasta)
  output: tuple path("p_${genome_fasta}"), path("m_${genome_fasta}")
  script:
  """
 #! /usr/bin/env bash
 #
 # === Inputs
 # genome_fasta = primary_assembly_merged.fasta
 #=== Outputs
 #  p_${genome_fasta}     # primary assembly
 #  m_${genome_fasta}     # mitochondrial assembly

  ${samtools_app} faidx ${genome_fasta}
  grep ">pri_" ${genome_fasta} | cut -f1 | sed 's/>//g' > pri.list
  ${samtools_app} faidx -r pri.list ${genome_fasta} > p_${genome_fasta}

  grep ">mit_" ${genome_fasta} | cut -f1 | sed 's/>//g' > mit.list
  ${samtools_app} faidx -r mit.list ${genome_fasta} > m_${genome_fasta}

  """

  stub:
  """
  touch p_${genome_fasta} m_${genome_fasta}
  """
}


process SPLIT_FILE {
  publishDir "${params.outdir}/${outdir}", mode:'copy'

  input:tuple val(outdir), path(genome_fasta)
  output: tuple path("p_${genome_fasta}"), path("a_${genome_fasta}"), path("m_${genome_fasta}")
  script:
  """
  #! /usr/bin/env bash

  # === Inputs
  # genome_fasta = primary_assembly_merged.fasta
  # === Outputs
  # p_${genome_fasta}     # primary assembly
  # a_${genome_fasta}     # alternative assembly
  # m_${genome_fasta}     # mitochondrial assembly

  ${samtools_app} faidx ${genome_fasta}
  grep ">pri_" ${genome_fasta} | cut -f1 | sed 's/>//g' > pri.list
  ${samtools_app} faidx -r pri.list ${genome_fasta} > p_${genome_fasta}

  grep ">mit_" ${genome_fasta} | cut -f1 | sed 's/>//g' > mit.list
  ${samtools_app} faidx -r mit.list ${genome_fasta} > m_${genome_fasta}

  grep ">alt_" ${genome_fasta} | cut -f1 | sed 's/>//g' > alt.list
  ${samtools_app} faidx -r alt.list ${genome_fasta} > a_${genome_fasta}
  """

  stub:
  """
  touch p_${genome_fasta} a_${genome_fasta} m_${genome_fasta}
  """
}

process RENAME_FILE {
  publishDir "${params.outdir}/00_Preprocess/", mode:'copy'

  input: tuple path(filename), val(newname)
  output: path("$newname")
  script:
  """
  #! /usr/bin/env bash
  mv ${filename} ${newname}
  """

  stub:
  """
  touch ${newname}
  """
}

// Used by both arrow and freebayes
process create_windows {
  publishDir "${params.outdir}/${outdir}", mode: 'symlink'
  input: tuple val(outdir), path(assembly_fasta)
  output: tuple path("*.fai"), path("win.txt")
  script:
  """
  #! /usr/bin/env bash
  ${samtools_app} faidx ${assembly_fasta}
  cat ${assembly_fasta}.fai | awk '{print \$1 ":0-" \$2}' > win.txt
  """

  stub:
  """
  touch ${assembly_fasta}.fai win.txt
  echo "1-10" >> win.txt
  echo "2-10" >> win.txt
  echo "3-10" >> win.txt
  """
}

process combineVCF {
  publishDir "${params.outdir}/${outdir}", mode: 'symlink'
  input: tuple val(outdir), path(vcfs)
  output: tuple val("$outdir"), path("*_consensus.vcf")
  script:
  """
#! /usr/bin/env bash

OUTNAME=`echo $outdir | sed 's:/:_:g'`
OUTVCF=\${OUTNAME}_consensus.vcf

# Merge by sections (1) headers up to contig, (2) all contig headers, (3) headers post contigs, (4) snp data
cat ${vcfs.get(0)} | sed -n '1,/##reference=/'p > \$OUTVCF
cat ${vcfs} | grep -h "##contig=" >> \$OUTVCF
cat ${vcfs.get(0)} | sed -n '/##INFO=/,/#CHROM/'p >> \$OUTVCF
cat ${vcfs} | grep -hv "#" >> \$OUTVCF
"""

  stub:
  """
  OUTNAME=`echo $outdir | sed 's:/:_:g'`
  touch \${OUTNAME}_consensus.vcf
  """
}

process meryl_genome {
  publishDir "${params.outdir}/${outdir}/merfin", mode: 'symlink'
  input: tuple val(outdir), val(k), path(illumina_read)
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

process merfin_polish {
  publishDir "${params.outdir}/${outdir}/merfin", mode: 'symlink'
  input: tuple val(outdir), path(vcf), path(genome_fasta), path(genome_meryl), val(peak), path(meryldb)
  output: tuple val("$outdir"), path("*merfin.polish.vcf")
  script:
  """
  #! /usr/bin/env bash

  OUTNAME=`echo $outdir | sed 's:/:_:g'`

  ${merfin_app} -polish \
    -sequence ${genome_fasta} \
    -seqmers ${genome_meryl} \
    -readmers ${meryldb} \
    -peak ${peak} \
    -vcf ${vcf} \
    -output \${OUTNAME}_merfin \
    ${merfin_params}

  """

  stub:
  """
  OUTNAME=`echo $outdir | sed 's:/:_:g'`
  touch \${OUTNAME}_merfin.polish.vcf
  """
}

process vcf_to_fasta {
  publishDir "${params.outdir}/${outdir}", mode: 'copy'
  input: tuple val(outdir), path(vcf), path(genome_fasta)
  output: path("*_consensus.fasta")
  script:
  """
  #! /usr/bin/env bash

  OUTNAME=`echo $outdir | sed 's:/:_:g'`

  ${bcftools_app} view -Oz ${vcf} > ${vcf}.gz
  ${bcftools_app} index ${vcf}.gz
  ${bcftools_app} consensus ${vcf}.gz -f ${genome_fasta} -Hla > \${OUTNAME}_consensus.fasta
  """

  stub:
  """
  OUTNAME=`echo $outdir | sed 's:/:_:g'`
  touch \${OUTNAME}_consensus.fasta
  """
}

process bam_to_fasta {
  publishDir "${params.outdir}/00_Preprocess", mode: 'copy'
  input: path(bam)
  output: path("${bam.simpleName}.fasta")
  script:
  """
  #! /usr/bin/env bash
  ${samtools_app} fasta ${bam} > ${bam.simpleName}.fasta
  """

  stub:
  """
  touch ${bam.simpleName}.fasta
  """

}
