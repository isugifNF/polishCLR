#! /usr/bin/env nextflow

nextflow.enable.dsl=2

def helpMessage() {
  log.info isuGIFHeader()
  log.info """
   Usage:
   The typical command for running the pipeline are as follows:
   nextflow run main.nf --primary_assembly "*fasta" --illumina_reads "*{1,2}.fastq.bz2" --pacbio_reads "*_subreads.bam" -resume

   Mandatory arguments:
   --primary_assembly             genome assembly fasta file to polish
   --illumina_reads               paired end illumina reads, to be used for Merqury QV scores, and freebayes polish primary assembly
   --pacbio_reads                 pacbio reads in bam format, to be used to arrow polish primary assembly

   Optional modifiers
   --k                            kmer to use in MerquryQV scoring [default:21]
   --same_specimen                if illumina and pacbio reads are from the same specimin [default: true]. 
   --falcon_unzip                 if primary assembly has already undergone falcon unzip [default: false]. If true, will Arrow polish once instead of twice.

   Optional arguments:
   --outdir                       Output directory to place final output [default: 'PolishCLR_Results']
   --clusterOptions               Cluster options for slurm or sge profiles [default slurm: '-N 1 -n 40 -t 04:00:00'; default sge: ' ']
   --threads                      Number of CPUs to use during each job [default: 40]
   --queueSize                    Maximum number of jobs to be queued [default: 50]
   --account                      Some HPCs require you supply an account name for tracking usage.  You can supply that here.
   --help                         This usage statement.
  """
}
//   --bzip                         if illumina paired reads are in bz2 format [default: true]. If true, will convert to gz.

// Show help message
if (params.help) {
  helpMessage()
  exit 0
}

if (!params.primary_assembly) {
  helpMessage()
  exit 0
}

if (!params.illumina_reads) {
  helpMessage()
  exit 0
}

if (!params.pacbio_reads) {
  helpMessage()
  exit 0
}

// 00 Preprocess: Unzip any bz2 files
process bz_to_gz {
    publishDir "${params.outdir}/00_Preprocess", mode: 'copy'
    input:tuple val(readname), path(illumina_reads)
    output: tuple val(readname), path("*.gz")
    script:
    """
    #! /usr/bin/env bash
    PROC=\$(((`nproc`-1)/2+1))
    ${parallel_app} -j 2 "${bzcat_app} {1} | ${pigz_app} -p \${PROC} > {1/.}.gz" ::: *.bz2
    """
}

// 01 Merqury: Quality value of primary assembly
process meryl_count {
    publishDir "${params.outdir}/01_MerquryQV", mode: 'symlink'
    input: tuple val(k), path(illumina_read)
    output: path("*.meryl")
    script:
    """
    #! /usr/bin/env bash
    $meryl_app count k=${k} output ${illumina_read.simpleName}.meryl ${illumina_read}
    """
}

process meryl_union {
    publishDir "${params.outdir}/01_MerquryQV", mode: 'copy'
    input: path(illumina_meryls)
    output: path("illumina.meryl")
    script:
    """
    #! /usr/bin/env bash
    $meryl_app union-sum output illumina.meryl ${illumina_meryls}
    """
}

workflow mk_MerylDB {
  take:
    k_ch
    illumina_ch
  main:
    meryldb_ch = k_ch | combine(illumina_ch) | meryl_count | collect | meryl_union
  emit:
    meryldb_ch    
}

process MerquryQV {
    publishDir "${params.outdir}/01_MerquryQV", mode: 'copy'
    input: tuple path(illumina_db), path(assembly_fasta)
    output: path("*")
    script:
    """
    #! /usr/bin/env bash
    $MERQURY/merqury.sh $illumina_db $assembly_fasta ${assembly_fasta.simpleName}
    """
}

workflow QV_01 {
  take:
    meryldb_ch
    fasta_ch
  main:
    qv_ch = meryldb_ch | combine(fasta_ch) | MerquryQV
  emit:
    qv_ch
}

// 1st Arrow Polish
process pbmm2_index {
    publishDir "${params.outdir}/02_ArrowPolish", mode: 'symlink'
    input:path(assembly_fasta)
    output:tuple path("$assembly_fasta"), path("*.mmi")
    script:
    """
    #! /usr/bin/env bash
    ${pbmm2_app} index ${assembly_fasta} ${assembly_fasta}.mmi
    """
}

process pbmm2_align {
    publishDir "${params.outdir}/02_ArrowPolish", mode: 'symlink'
    input:tuple path(assembly_fasta), path(assembly_mmi), path(pacbio_read)
    output: tuple path("*.bam"), path("*.bai")
    script:
    """
    #! /usr/bin/env bash
    PROC=\$((`nproc` /2+1))
    mkdir tmp
    ${pbmm2_app} align -j \$PROC ${assembly_fasta} ${pacbio_read} | ${samtools_app} sort -T tmp -m 8G --threads 8 - > ${pacbio_read.simpleName}_aln.bam
    ${samtools_app} index -@ \${PROC} ${pacbio_read.simpleName}_aln.bam
    """
}

process create_windows {
    publishDir "${params.outdir}/02_ArrowPolish", mode: 'symlink'
    input: path(assembly_fasta)
    output: tuple path("*.fai"), path("win_02.txt")
    shell:
    """
    #! /usr/bin/env bash
    ${samtools_app} faidx ${assembly_fasta}
    cat ${assembly_fasta}.fai | awk '{print \$1 ":0-" \$2}' > win.txt
    """
}

process gcc_Arrow {
    publishDir "${params.outdir}/02_ArrowPolish", mode: 'symlink'
    input: tuple val(window), path(assembly_fasta), path(assembly_fai), path(pacbio_bam), path(pacbio_bai)
    output: tuple path("*.fasta"), path("*.vcf")
    script:
    """
    #! /usr/bin/env bash
    PROC=\$((`nproc` /2+1))
    ${gcpp_app} --algorithm=arrow \
      -x 10 -X 120 -q 0 \
      -j \${PROC} -w \"$window\" \
      -r ${assembly_fasta} ${pacbio_bam} \
      -o ${window.replace(':','_').replace('|','_')}.vcf,${window.replace(':','_').replace('|','_')}.fasta
    """
}

// add nproc
process merge_consensus {
    publishDir "${params.outdir}/02_ArrowPolish", mode: 'symlink'
    input: path(windows_fasta)
    output: path("consensus.fasta")
    script:
    """
    #! /usr/bin/env bash
    cat ${windows_fasta} > consensus.fasta
    """
}

// 2nd Merqury QV value
workflow Arrow_02 {
  take:
    asm_ch
    pac_ch
  main:
    // Step 2: Arrow Polish with PacBio reads
    asm_ch | pbmm2_index | combine(pac_ch) | pbmm2_align
    fai_ch = asm_ch | create_windows | map { n -> n.get(0) }
    create_windows.out | 
      map { n -> n.get(1) } | 
      splitText() {it.trim()} |
      combine(asm_ch) | combine(fai_ch) | combine(pbmm2_align.out) | 
      gcc_Arrow |
      map { n -> n.get(0)} |
      collect |
      merge_consensus
    new_asm_ch = merge_consensus.out   // <= New Assembly
  emit:
    new_asm_ch
}

workflow QV_03 {
  take:
    meryldb_ch
    fasta_ch
  main:
    qv_ch = meryldb_ch | combine(fasta_ch) | MerquryQV
  emit:
    qv_ch
}

// 2nd Arrow Polish (skip if falcon unzip)
workflow Arrow_04 {
  take:
    asm_ch
    pac_ch
  main:
    // Step 2: Arrow Polish with PacBio reads
    asm_ch | pbmm2_index | combine(pac_ch) | pbmm2_align
    fai_ch = asm_ch | create_windows | map { n -> n.get(0) }
    create_windows.out | 
      map { n -> n.get(1) } | 
      splitText() {it.trim()} |
      combine(asm_ch) | combine(fai_ch) | combine(pbmm2_align.out) | 
      gcc_Arrow |
      map { n -> n.get(0)} |
      collect |
      merge_consensus
    new_asm_ch = merge_consensus.out   // <= New Assembly
  emit:
    new_asm_ch
}

workflow QV_05 {
  take:
    meryldb_ch
    fasta_ch
  main:
    qv_ch = meryldb_ch | combine(fasta_ch) | MerquryQV
  emit:
    qv_ch
}

// 1st FreeBayes Polish
process align_shortreads {
    publishDir "${params.outdir}/04_FreeBayesPolish", mode: 'symlink'
    input: tuple path(assembly_fasta), path(illumina_one), path(illumina_two)
    output: tuple path("*.bam"), path("*.bai")
    script:
    """
    #! /usr/bin/env bash
    PROC=\$((`nproc` /2+1))
    mkdir tmp
    ${bwamem2_app} index ${assembly_fasta}
    ${bwamem2_app} mem -SP -t \${PROC} ${assembly_fasta} ${illumina_one} ${illumina_two} |
      ${samtools_app} sort -T tmp -m 8G --threads 4 - > ${illumina_one.simpleName}_aln.bam
    ${samtools_app} index -@ \${PROC} ${illumina_one.simpleName}_aln.bam
    """
}
// -m 5G -@ 36

process freebayes {
    publishDir "${params.outdir}/04_FreeBayesPolish", mode: 'symlink'
    input: tuple val(window), path(assembly_fasta), path(assembly_fai), path(illumina_bam), path(illumina_bai)
    output: path("*.vcf")
    script:
    """
    #! /usr/bin/env bash
    ${freebayes_app} \
      --region \"${window}\" \
      --min-mapping-quality 0 \
      --min-coverage 3 \
      --min-supporting-allele-qsum 0 \
      --ploidy 2 \
      --min-alternate-fraction 0.2 \
      --max-complex-gap 0 \
      --bam ${illumina_bam} \
      --vcf ${illumina_bam.simpleName}_${window.replace(':','_').replace('|','_')}.vcf \
      --fasta-reference ${assembly_fasta}
    """
}

process combineVCF {
    publishDir "${params.outdir}/04_FreeBayesPolish", mode: 'symlink'
    input: path(vcfs)
    output: path("consensus_04.vcf")
    script:
    """
    #! /usr/bin/env bash
    cat ${vcfs.get(0)} | grep "^#" > consensus_04.vcf
    cat ${vcfs} | grep -v "^#" >> consensus_04.vcf
    """
}

process vcf_to_fasta {
    publishDir "${params.outdir}/04_FreeBayesPolish", mode: 'symlink'
    input: tuple path(vcf), path(genome_fasta)
    output: path("consensus_04.fasta")
    script:
    """
    #! /usr/bin/env bash
    ${bcftools_app} view -Oz ${vcf} > ${vcf.simpleName}.gz
    ${bcftools_app} index ${vcf.simpleName}.gz
    ${bcftools_app} consensus ${vcf.simpleName}.gz -f ${genome_fasta} -H 1 > consensus_04.fasta
    """
}

workflow Freebayes_06 {
  take:
    assembly_ch
    illumina_ch
  main:
    assembly_ch | combine(illumina_ch.collect()) | align_shortreads 
    fai_ch = assembly_ch | create_windows | map { n -> n.get(0) }
    create_windows.out | 
      map { n -> n.get(1) } | 
      splitText() {it.trim()} |
      combine(assembly_ch) | combine(fai_ch) | combine(align_shortreads.out) |
      freebayes |
      collect |
      combineVCF | 
      combine(assembly_ch) |
      vcf_to_fasta
    new_asm_ch = vcf_to_fasta.out         // <= New assembly
  emit:
    new_asm_ch
}


// 3rd Merqury QV value
process MerquryQV_05 {
    publishDir "${params.outdir}/05_MerquryQV", mode: 'symlink'
    input: tuple path(illumina_db), path(assembly_fasta)
    output: path("*")
    script:
    """
    #! /usr/bin/env bash
    $MERQURY/merqury.sh $illumina_db $assembly_fasta ${assembly_fasta.simpleName}
    """
}

// 2nd Freebayes polish
process align_shortreads_06 {
    publishDir "${params.outdir}/06_FreeBayesPolish", mode: 'symlink'
    input: tuple path(assembly_fasta), path(illumina_one), path(illumina_two)
    output: tuple path("*.bam"), path("*.bai")
    script:
    """
    #! /usr/bin/env bash
    PROC=\$((`nproc` /2+1))
    mkdir tmp
    bwa-mem2 index ${assembly_fasta}
    bwa-mem2 mem -SP -t \${PROC} ${assembly_fasta} ${illumina_one} ${illumina_two} |
      samtools sort -T tmp -m 8G --threads 4 - > ${illumina_one.simpleName}_06_aln.bam
    samtools index -@ \${PROC} ${illumina_one.simpleName}_06_aln.bam
    """
}
// -m 5G -@ 36

process create_windows_06 {
    publishDir "${params.outdir}/06_FreeBayesPolish", mode: 'symlink'
    input: path(assembly_fasta)
    output: tuple path("*.fai"), path("win_06.txt")
    shell:
    """
    #! /usr/bin/env bash
    samtools faidx ${assembly_fasta}
    cat ${assembly_fasta}.fai | awk '{print \$1 ":0-" \$2}' > win_06.txt
    """
}

process freebayes_06 {
    publishDir "${params.outdir}/06_FreeBayesPolish", mode: 'symlink'
    input: tuple val(window), path(assembly_fasta), path(assembly_fai), path(illumina_bam), path(illumina_bai)
    output: path("*.vcf")
    script:
    """
    #! /usr/bin/env bash
    freebayes \
      --region \"${window}\" \
      --min-mapping-quality 0 \
      --min-coverage 3 \
      --min-supporting-allele-qsum 0 \
      --ploidy 2 \
      --min-alternate-fraction 0.2 \
      --max-complex-gap 0 \
      --bam ${illumina_bam} \
      --vcf ${illumina_bam.simpleName}_${window.replace(':','_').replace('|','_')}_06.vcf \
      --fasta-reference ${assembly_fasta}
    """
}

process combineVCF_06 {
    publishDir "${params.outdir}/06_FreeBayesPolish", mode: 'symlink'
    input: path(vcfs)
    output: path("consensus_06.vcf")
    script:
    """
    #! /usr/bin/env bash
    cat ${vcfs.get(0)} | grep "^#" > consensus_06.vcf
    cat ${vcfs} | grep -v "^#" >> consensus_06.vcf
    """
}

process vcf_to_fasta_06 {
    publishDir "${params.outdir}/06_FreeBayesPolish", mode: 'symlink'
    input: tuple path(vcf), path(genome_fasta)
    output: path("final_polished_assembly.fasta")
    script:
    """
    #! /usr/bin/env bash
    bcftools view -Oz ${vcf} > ${vcf}.gz
    bcftools index ${vcf.simpleName}.vcf.gz
    bcftools consensus ${vcf.simpleName}.vcf.gz -f ${genome_fasta} -H 1 > final_polished_assembly.fasta
    """
}

// 3rd Merqury QV value
process MerquryQV_07 {
    publishDir "${params.outdir}/07_MerquryQV", mode: 'symlink'
    input: tuple path(illumina_db), path(assembly_fasta)
    output: path("*")
    script:
    """
    #! /usr/bin/env bash
    $MERQURY/merqury.sh $illumina_db $assembly_fasta ${assembly_fasta.simpleName}
    """
}

workflow {
    // Setup input channels, starting assembly (asm), Illumina reads (ill), and pacbio reads (pac)
    asm_ch = channel.fromPath(params.primary_assembly, checkIfExists:true)
    ill_ch = channel.fromFilePairs(params.illumina_reads, checkIfExists:true) 
    pac_ch = channel.fromPath(params.pacbio_reads, checkIfExists:true)
    k_ch   = channel.of(params.k)

    // Step 0: Preprocess illumina files from bz2 to gz files
//    if(params.bzip2){     // Instead of a flag, auto detect, however it must be in the pattern, * will fail
    if(params.illumina_reads =~ /bz2$/){
      pill_ch = ill_ch | bz_to_gz | map { n -> n.get(1) } | flatten
    }else{
      pill_ch = ill_ch | map { n -> n.get(1) } | flatten
    }

    // Step 1: Check quality of assembly with Merqury
    meryldb_ch = mk_MerylDB(k_ch, pill_ch)
    QV_01(meryldb_ch, asm_ch)

    // Step 2: Arrow Polish with PacBio reads
    asm_arrow_ch = Arrow_02(asm_ch, pac_ch)

    // Step 3: Check quality of new assembly with Merqury (turns out we can reuse the illumina database)
    QV_03(meryldb_ch, asm_arrow_ch)

    // if the primary assembly came from falcon unzip, skip the 2nd arrow polish
    if(!params.falcon_unzip) {
      // Step 2b: Arrow Polish with PacBio reads
      asm_arrow2_ch = Arrow_04(asm_arrow_ch, pac_ch)
      // Step 3b: Check quality of new assembly with Merqury (turns out we can reuse the illumina database)
      QV_05(meryldb_ch, asm_arrow2_ch)
    } else {
      asm_arrow2_ch = asm_arrow_ch
    }

    // Step 4: FreeBayes Polish with Illumina reads

    asm2b_ch | combine(pill_ch.collect()) | align_shortreads_04 | view
    fai2_ch = asm2b_ch | create_windows_04 | map { n -> n.get(0) } | view
    create_windows_04.out | 
      map { n -> n.get(1) } | 
      splitText() {it.trim()} |
      combine(asm2b_ch) | combine(fai2_ch) | combine(align_shortreads_04.out) |
      freebayes_04 |
      collect |
      combineVCF_04 | 
      combine(asm2b_ch) |
      vcf_to_fasta_04
    asm3_ch = vcf_to_fasta_04.out         // <= New assembly
  
    // Step 5: Check quality of assembly with Merqury
    meryl_union_01.out | combine(asm3_ch) | MerquryQV_05

    // Step 6: FreeBayes Polish 2nd time
    asm3_ch | combine(pill_ch.collect()) | align_shortreads_06 | view
    fai3_ch = asm3_ch | create_windows_06 | map { n -> n.get(0) } | view
    create_windows_06.out | 
      map { n -> n.get(1) } | 
      splitText() {it.trim()} |
      combine(asm3_ch) | combine(fai3_ch) | combine(align_shortreads_06.out) |
      freebayes_06 |
      collect |
      combineVCF_06 | 
      combine(asm3_ch) |
      vcf_to_fasta_06
    asm4_ch = vcf_to_fasta_06.out         // <= New assembly

    // Step 7: Check quality of assembly with Merqury
    meryl_union_01.out | combine(asm4_ch) | MerquryQV_07
}

def isuGIFHeader() {
  // Log colors ANSI codes
  c_reset = params.monochrome_logs ? '' : "\033[0m";
  c_dim = params.monochrome_logs ? '' : "\033[2m";
  c_black = params.monochrome_logs ? '' : "\033[1;90m";
  c_green = params.monochrome_logs ? '' : "\033[1;92m";
  c_yellow = params.monochrome_logs ? '' : "\033[1;93m";
  c_blue = params.monochrome_logs ? '' : "\033[1;94m";
  c_purple = params.monochrome_logs ? '' : "\033[1;95m";
  c_cyan = params.monochrome_logs ? '' : "\033[1;96m";
  c_white = params.monochrome_logs ? '' : "\033[1;97m";
  c_red = params.monochrome_logs ? '' :  "\033[1;91m";

  return """    -${c_dim}--------------------------------------------------${c_reset}-
  ${c_white}                                ${c_red   }\\\\------${c_yellow}---//       ${c_reset}
  ${c_white}  ___  ___        _   ___  ___  ${c_red   }  \\\\---${c_yellow}--//        ${c_reset}
  ${c_white}   |  (___  |  | / _   |   |_   ${c_red   }    \\-${c_yellow}//         ${c_reset}
  ${c_white}  _|_  ___) |__| \\_/  _|_  |    ${c_red  }    ${c_yellow}//${c_red  } \\        ${c_reset}
  ${c_white}                                ${c_red   }  ${c_yellow}//---${c_red  }--\\\\       ${c_reset}
  ${c_white}                                ${c_red   }${c_yellow}//------${c_red  }---\\\\       ${c_reset}
  ${c_cyan}  isugifNF/polishCLR  v${workflow.manifest.version}       ${c_reset}
  -${c_dim}--------------------------------------------------${c_reset}-
  """.stripIndent()
}
