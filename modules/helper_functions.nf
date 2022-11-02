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

process PREFIX_FASTA {
  input: tuple path(fasta), val(prefix)
  output: path("${prefix}_$fasta")
  script:
  """
  #! /usr/bin/env bash
  cat ${fasta} | sed 's/>/>${prefix}_/g' > ${prefix}_${fasta}
  """
  stub:
  """
  touch ${prefix}_${fasta}
  """
}

process CONCATINATE_FASTA {
  input: tuple path(fastas), val(concatinated_fasta)
  output: path("${concatinated_fasta}")
  script:
  """
  #! /usr/bin/env bash

  echo "Expecting output: ${concatinated_fasta}\n"

  touch temp.fasta
  for FILE in "${fastas}"; do
    cat \${FILE} >> temp.fasta
    echo "" >> temp.fasta
  done
  grep -v "^\$" temp.fasta > ${concatinated_fasta}
  
  ls -ltr
  """
  stub:
  """
  touch "${concatinated_fasta}"
  """
}

process SPLIT_FASTA {
  input: path(fasta)
  output: path("*_${fasta}")
  script:
  """
  #! /usr/bin/env bash
  ${samtools_app} faidx ${fasta}
  for PREFIX in "pri" "alt" "mit" ; do
    grep "\${PREFIX}_" ${fasta}.fai | cut -f1  > \${PREFIX}.list
    if [[ -s "\${PREFIX}.list" ]]; then
      ${samtools_app} faidx -r \${PREFIX}.list ${fasta} > \${PREFIX}_${fasta}
    fi
  done
  """
  stub:
  """
  touch pri_${fasta} alt_${fasta} mit_${fasta}
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
