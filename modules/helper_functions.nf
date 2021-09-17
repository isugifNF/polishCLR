#! /usr/bin/env nextflow

nextflow.enable.dsl=2

// 00 Preprocess: Unzip any bz2 files
process bz_to_gz {
  publishDir "${params.outdir}/00_Preprocess", mode: 'copy'
  input:tuple val(readname), path(illumina_reads)
  output: tuple val(readname), path("*.gz")
  script:
  template 'bz_to_gz.sh'

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
  template 'merge_file.sh'

  stub:
  """
  touch ${primary_assembly.simpleName}_merged.fasta
  """

}

process MERGE_FILE_TRIO {
  publishDir "${params.outdir}/00_Preprocess", mode: 'copy'
  input: tuple path(primary_assembly), path(mito_assembly)
  output: path("${primary_assembly.simpleName}_merged.fasta")
  script:
  template 'merge_file_trio.sh'

  stub:
  """
  touch ${primary_assembly.simpleName}_merged.fasta
  """
}

process SPLIT_FILE_03p {
  input: path(genome_fasta)
  output: tuple path("pat_${genome_fasta}"), path("mit_${genome_fasta}")
  script:
  template 'split_file_pat.sh'

  stub:
  """
  touch pat_${genome_fasta} mit_${genome_fasta}
  """
}

process SPLIT_FILE_03m {
  input: path(genome_fasta)
  output: tuple path("mat_${genome_fasta}"), path("mit_${genome_fasta}")
  script:
  template 'split_file_mat.sh'

  stub:
  """
  touch mat_${genome_fasta} mit_${genome_fasta}
  """
}

process SPLIT_FILE {
  input: path(genome_fasta)
  output: tuple path("p_${genome_fasta}"), path("a_${genome_fasta}"), path("m_${genome_fasta}")
  script:
  template 'split_file.sh'

  stub:
  """
  touch p_${genome_fasta} a_${genome_fasta} m_${genome_fasta}
  """
}


// Used by both arrow and freebayes
process create_windows {
  publishDir "${params.outdir}/${outdir}", mode: 'symlink'
  input: tuple val(outdir), path(assembly_fasta)
  output: tuple path("*.fai"), path("win.txt")
  script:
  template 'create_windows.sh'

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
  output: tuple val("$outdir"), path("${outdir}_consensus.vcf")
  script:
  template 'combineVCF.sh'

  stub:
  """
  touch ${outdir}_consensus.vcf
  """
}

process meryl_genome {
  publishDir "${params.outdir}/${outdir}/merfin", mode: 'symlink'
  input: tuple val(outdir), val(k), path(illumina_read)
  output: path("*.meryl")
  script:
  template 'meryl_count.sh'

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
  template 'merfin_polish.sh'

  stub:
  """
  touch ${outdir}_merfin.polish.vcf
  """
}

process vcf_to_fasta {
  publishDir "${params.outdir}/${outdir}", mode: 'copy'
  input: tuple val(outdir), path(vcf), path(genome_fasta)
  output: path("${outdir}_consensus.fasta")
  script:
  template 'vcf_to_fasta.sh'

  stub:
  """
  touch ${outdir}_consensus.fasta
  """
}