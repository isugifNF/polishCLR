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

   Optional configuration arguments
   --parallel_app                 Link to parallel executable [default: 'parallel']
   --bzcat_app                    Link to bzcat executable [default: 'bzcat']
   --pigz_app                     Link to pigz executable [default: 'pigz']
   --meryl_app                    Link to meryl executable [default: 'meryl']
   --merqury_sh                   Link to merqury script [default: '\$MERQURY/merqury.sh']
   --pbmm2_app                    Link to pbmm2 executable [default: 'pbmm2']
   --samtools_app                 Link to samtools executable [default: 'samtools']
   --gcpp_app                     Link to gcpp executable [default: 'gcpp']
   --bwamem2_app                  Link to bwamem2 executable [default: 'bwa-mem2']
   --freebayes_app                Link to freebayes executable [default: 'freebayes']
   --bcftools_app                 Link to bcftools executable [default: 'bcftools']

   Optional arguments:
   --outdir                       Output directory to place final output [default: 'PolishCLR_Results']
   --clusterOptions               Cluster options for slurm or sge profiles [default slurm: '-N 1 -n 40 -t 04:00:00'; default sge: ' ']
   --threads                      Number of CPUs to use during each job [default: 40]
   --queueSize                    Maximum number of jobs to be queued [default: 50]
   --account                      Some HPCs require you supply an account name for tracking usage.  You can supply that here.
   --help                         This usage statement.
  """
}

// Show help message
if (params.help || !params.primary_assembly || !params.illumina_reads || !params.pacbio_reads ) {
  helpMessage()
  exit 0
}

// 00 Preprocess: Unzip any bz2 files
process bz_to_gz {
    publishDir "${params.outdir}/00_Preprocess", mode: 'copy'
    input:tuple val(readname), path(illumina_reads)
    output: tuple val(readname), path("*.gz")
    script:
    template 'bz_to_gz.sh'
}

// 01 Merqury: Quality value of primary assembly
process meryl_count {
    publishDir "${params.outdir}/01_MerquryQV", mode: 'symlink'
    input: tuple val(k), path(illumina_read)
    output: path("*.meryl")
    script:
    template 'meryl_count.sh'
}

process meryl_union {
    publishDir "${params.outdir}/01_MerquryQV", mode: 'copy'
    input: path(illumina_meryls)
    output: path("illumina.meryl")
    script:
    template 'meryl_union.sh'
}

process MerquryQV_01 {
    publishDir "${params.outdir}/01_MerquryQV", mode: 'copy'
    input: tuple path(illumina_db), path(assembly_fasta)
    output: path("*")
    script:
    template 'merquryqv.sh'
}

// 1st Arrow Polish
process create_windows {
    publishDir "${params.outdir}/0${i}_ArrowPolish", mode: 'symlink'
    input: tuple val(i), path(assembly_fasta)
    output: tuple path("*.fai"), path("win.txt")
    script:
    template 'create_windows.sh'
}

process pbmm2_index {
    publishDir "${params.outdir}/0${i}_ArrowPolish", mode: 'symlink'
    input: tuple val(i), path(assembly_fasta)
    output: tuple val("$i"), path("$assembly_fasta"), path("*.mmi")
    script:
    template 'pbmm2_index.sh'
}

process pbmm2_align {
    publishDir "${params.outdir}/0${i}_ArrowPolish", mode: 'symlink'
    input:tuple val(i), path(assembly_fasta), path(assembly_mmi), path(pacbio_read)
    output: tuple val("$i"), path("*.bam"), path("*.bai")
    script:
    template 'pbmm2_align.sh'
}

process gcc_Arrow {
    publishDir "${params.outdir}/0${i}_ArrowPolish", mode: 'symlink'
    input: tuple val(i), path(pacbio_bam), path(pacbio_bai),  path(assembly_fasta), path(assembly_fai), val(window)
    output: tuple val("$i"), path("*.fasta"), path("*.vcf")
    script:
    template 'gcc_arrow.sh'
}

process merge_consensus {
    publishDir "${params.outdir}/0${i}_ArrowPolish", mode: 'copy'
    input: tuple val(i), path(windows_fasta)
    output: path("${i}_consensus.fasta")
    script:
    template 'merge_consensus.sh'
}

workflow ARROW_02 {
  take:
    asm_ch
    pac_ch
  main:
    win_ch = channel.of("2") | combine(asm_ch) | create_windows | 
      map { n -> n.get(1) } | splitText() {it.trim() }
    fai_ch = create_windows.out | map { n -> n.get(0) }

    newasm_ch = channel.of("2") | combine(asm_ch) | pbmm2_index | combine(pac_ch) | pbmm2_align |
      combine(asm_ch) | combine(fai_ch) | combine(win_ch) | gcc_Arrow | 
      map { n -> [ n.get(0), n.get(1) ] } | groupTuple | merge_consensus
  
  emit:
    newasm_ch
}

process MerquryQV_03 {
    publishDir "${params.outdir}/03_MerquryQV", mode: 'copy'
    input: tuple path(illumina_db), path(assembly_fasta)
    output: path("*")
    script:
    template 'merquryqv.sh'
}

// 2nd Arrow run with merfin
process reshape_arrow {
    publishDir "${params.outdir}/0${i}_ArrowPolish", mode: 'copy'
    input: tuple val(i), path(vcf), path(asm)
    output: path("${i}_merged.reshaped.vcf.gz")
    script:
    template 'reshape_arrow.sh'
}

workflow ARROW_04 {
  take:
    asm_ch
    pac_ch
  main:
    win_ch = channel.of("4") | combine(asm_ch) | create_windows | 
      map { n -> n.get(1) } | splitText() {it.trim() }
    fai_ch = create_windows.out | map { n -> n.get(0) }

    newasm_ch = channel.of("4") | combine(asm_ch) | pbmm2_index | combine(pac_ch) | pbmm2_align |
      combine(asm_ch) | combine(fai_ch) | combine(win_ch) | gcc_Arrow | 
      map { n -> [ n.get(0), n.get(2) ] } | groupTuple |
      combine(asm_ch) |  reshape_arrow  //| 
//      groupTuple | merge_consensus
  
  emit:
    newasm_ch
}

process MerquryQV_05 {
    publishDir "${params.outdir}/05_MerquryQV", mode: 'copy'
    input: tuple path(illumina_db), path(assembly_fasta)
    output: path("*")
    script:
    template 'merquryqv.sh'
}

// 1st FreeBayes Polish
process align_shortreads {
    publishDir "${params.outdir}/0${i}_FreeBayesPolish", mode: 'symlink'
    input: tuple val(i), path(assembly_fasta), path(illumina_one), path(illumina_two)
    output: tuple val("$i"), path("*.bam"), path("*.bai")
    script:
    template 'align_shortreads.sh'
}
// // -m 5G -@ 36
//
process freebayes {
    publishDir "${params.outdir}/0${i}_FreeBayesPolish", mode: 'symlink'
    input: tuple val(i), path(illumina_bam), path(illumina_bai), path(assembly_fasta), path(assembly_fai), val(window)
    output: tuple val(i), path("*.vcf")
    script:
    template 'freebayes.sh'
}

process combineVCF {
    publishDir "${params.outdir}/0${i}_FreeBayesPolish", mode: 'symlink'
    input: tuple val(i), path(vcfs)
    output: tuple val(i), path("${i}_consensus.vcf")
    script:
    template 'combineVCF.sh'
}

process jellyfish_peak {
    publishDir "${params.outdir}/0${i}_MerfinPolish", mode: 'copy'
    input: tuple val(i), path(illumina_reads)
    output: path("*")
    script:
    template 'jellyfish_peak.sh'
}

process merfin_polish {
    publishDir "${params.outdir}/0${i}_MerfinPolish", mode: 'symlink'
    input: tuple val(i), path(vcf), path(genome_fasta), val(peak), path(meryldb)
    output: path("*")
    script:
    template 'merfin_polish.sh'
}

process vcf_to_fasta {
    publishDir "${params.outdir}/0${i}_FreeBayesPolish", mode: 'symlink'
    input: tuple val(i), path(vcf), path(genome_fasta)
    output: path("${i}_consensus.fasta")
    script:
    template 'vcf_to_fasta.sh'
}

// First freebayes
workflow FREEBAYES_06 {
  take:
    asm_ch
    ill_ch
    merylDB_ch
  main:
    win_ch = channel.of("6") | combine(asm_ch) | create_windows | 
      map { n -> n.get(1) } | splitText() {it.trim() }
    fai_ch = create_windows.out | map { n -> n.get(0) }

    peak_ch = channel.of("6") | combine( ill_ch.collect() | map { n-> [n]} ) | jellyfish_peak | splitText() { it.trim()}

    new_asm_ch = channel.of("6") | combine(asm_ch) | combine(ill_ch.collect()) | align_shortreads |
      combine(asm_ch) | combine(fai_ch) | combine(win_ch) |
      freebayes | groupTuple |
      combineVCF |
      combine(asm_ch) | combine(peak_ch) | combine(merylDB_ch) //|
//      merfin_polish
//      vcf_to_fasta
  emit:
    new_asm_ch
}

process MerquryQV_07 {
    publishDir "${params.outdir}/07_MerquryQV", mode: 'copy'
    input: tuple path(illumina_db), path(assembly_fasta)
    output: path("*")
    script:
    template 'merquryqv.sh'
}

// 2nd Freebayes polish

workflow FREEBAYES_08 {
  take:
    asm_ch
    ill_ch
  main:
    win_ch = channel.of("8") | combine(asm_ch) | create_windows | 
      map { n -> n.get(1) } | splitText() {it.trim() }
    fai_ch = create_windows.out | map { n -> n.get(0) }

    new_asm_ch = channel.of("8") | combine(asm_ch) | combine(ill_ch.collect()) | align_shortreads |
      combine(asm_ch) | combine(fai_ch) | combine(win_ch) |
      freebayes | groupTuple |
      combineVCF |
      combine(asm_ch) |
      vcf_to_fasta
  emit:
    new_asm_ch
}

process MerquryQV_09 {
    publishDir "${params.outdir}/09_MerquryQV", mode: 'copy'
    input: tuple path(illumina_db), path(assembly_fasta)
    output: path("*")
    script:
    template 'merquryqv.sh'
}

workflow {
    // Setup input channels, starting assembly (asm), Illumina reads (ill), and pacbio reads (pac)
    asm_ch = channel.fromPath(params.primary_assembly, checkIfExists:true)
    ill_ch = channel.fromFilePairs(params.illumina_reads, checkIfExists:true)
    pac_ch = channel.fromPath(params.pacbio_reads, checkIfExists:true)
    k_ch   = channel.of(params.k) // Either passed in or autodetect (there's a script for this)

    // Step 0: Preprocess illumina files from bz2 to gz files
    // Instead of a flag, auto detect, however it must be in the pattern, * will fail
    if(params.illumina_reads =~ /bz2$/){
      pill_ch = ill_ch | bz_to_gz | map { n -> n.get(1) } | flatten
    }else{
      pill_ch = ill_ch | map { n -> n.get(1) } | flatten
    }

    // Step 1: Check quality of assembly with Merqury
    merylDB_ch = k_ch | combine(pill_ch) | meryl_count | collect | meryl_union 
    merylDB_ch | combine(asm_ch) | MerquryQV_01
    
    // Step 2: Arrow Polish with PacBio reads
    asm_arrow_ch = ARROW_02(asm_ch, pac_ch)

    // Step 3: Check quality of new assembly with Merqury 
    merylDB_ch | combine(asm_arrow_ch) | MerquryQV_03

    // if the primary assembly came from falcon unzip, skip the 2nd arrow polish
    if(!params.falcon_unzip) {
      // Step 4: Arrow Polish with PacBio reads
      asm_arrow2_ch = ARROW_04(asm_arrow_ch, pac_ch)
      // Step 5: Check quality of new assembly with Merqury 
      merylDB_ch | combine(asm_arrow2_ch) | MerquryQV_05
    } else {
      asm_arrow2_ch = asm_arrow_ch
    }
    asm_arrow2_ch | view
    
    // Step 6: FreeBayes Polish with Illumina reads
    asm_freebayes_ch = FREEBAYES_06(asm_arrow_ch, pill_ch, merylDB_ch)
/*    merylDB_ch | combine(asm_freebayes_ch) | MerquryQV_07

    // Step 8: FreeBayes Polish with Illumina reads
    asm_freebayes2_ch = FREEBAYES_08(asm_freebayes_ch, pill_ch)
    merylDB_ch | combine(asm_freebayes2_ch) | MerquryQV_09
*/
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
