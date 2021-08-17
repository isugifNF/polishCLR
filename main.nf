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
   --species                      if a string is given, rename the final assembly by species name [default:false]

   Optional modifiers
   --mito_assembly                if mitocondrial assembly file is provided, will be concatted to the primary assembly before polishing [default: false]
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
   --merfin_app                   Link to merfin executable [default: 'merfin']

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

// If user uses --profile, exit early. The param should be -profile (one hyphen)
if ( params.profile ) {
  helpMessage()
  println("Instead of --profile, use -profile")
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
    publishDir "${params.outdir}/01_QV", mode: 'symlink'
    input: tuple val(k), path(illumina_read)
    output: path("*.meryl")
    script:
    template 'meryl_count.sh'
}

process meryl_union {
    publishDir "${params.outdir}/01_QV", mode: 'copy'
    input: path(illumina_meryls)
    output: path("illumina.meryl")
    script:
    template 'meryl_union.sh'
}

process MerquryQV_01 {
    publishDir "${params.outdir}/01_QV/MerquryQV", mode: 'copy'
    publishDir "${params.outdir}/01_QV", mode: 'copy', pattern: "merqury.qv"
    input: tuple path(illumina_db), path(assembly_fasta)
    output: path("*")
    script:
    template 'merquryqv.sh'
}

// 01 bbstat: Length distribtions of initial assembly
process bbstat_01 {
    publishDir "${params.outdir}/01_QV/bbstat", mode: 'copy'
    input:  path(assembly_fasta)
    output: path("*")
    script:
    template 'bbstats.sh'
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
    publishDir "${params.outdir}/0${i}_ArrowPolish/gccruns", mode: 'symlink'
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
      map { n -> [ n.get(0), n.get(1) ] } | groupTuple | 
      merge_consensus
  
  emit:
    newasm_ch
}

process MerquryQV_03 {
    publishDir "${params.outdir}/03_QV/MerquryQV", mode: 'copy'
    publishDir "${params.outdir}/03_QV", mode: 'copy', pattern: "merqury.qv"
    input: tuple path(illumina_db), path(assembly_fasta)
    output: path("*")
    script:
    template 'merquryqv.sh'
}

process bbstat_03 {
    publishDir "${params.outdir}/03_QV/bbstat", mode: 'copy'
    input:  path(assembly_fasta)
    output: path("*")
    script:
    template 'bbstats.sh'
}

// 2nd Arrow run with merfin
process meryl_peak {
    publishDir "${params.outdir}/0${i}_ArrowPolish/merfin", mode: 'symlink'
    input: tuple val(i), path(illumina_meryl)
    output: path("peak.txt")
    script:
    template 'meryl_peak.sh'
}

process meryl_genome {
    publishDir "${params.outdir}/04_ArrowPolish/merfin", mode: 'symlink'
    input: tuple val(k), path(illumina_read)
    output: path("*.meryl")
    script:
    template 'meryl_count.sh'
}

process combineVCF_arrow {
    publishDir "${params.outdir}/0${i}_ArrowPolish", mode: 'symlink'
    input: tuple val(i), path(vcfs)
    output: tuple val(i), path("${i}_consensus.vcf")
    script:
    template 'combineVCF.sh'
}

process reshape_arrow {
    publishDir "${params.outdir}/0${i}_ArrowPolish/merfin", mode: 'symlink'
    input: tuple val(i), path(vcf)
    output: tuple val(i), path("*.reshaped.vcf.gz")
    script:
    template 'reshape_arrow.sh'
}

process merfin_polish_arrow {
    publishDir "${params.outdir}/0${i}_ArrowPolish/merfin", mode: 'symlink'
    input: tuple val(i), path(vcf), path(genome_fasta), path(genome_meryl), val(peak), path(meryldb)
    output: tuple val("$i"), path("*merfin.polish.vcf")
    script:
    template 'merfin_polish_arrow.sh'
}

process vcf_to_fasta_arrow {
    publishDir "${params.outdir}/0${i}_ArrowPolish", mode: 'symlink'
    input: tuple val(i), path(vcf), path(genome_fasta)
    output: path("${i}_consensus.fasta")
    script:
    template 'vcf_to_fasta.sh'
}

workflow ARROW_04 {
  take:
    asm_ch
    pac_ch
    merylDB_ch
    
  main:
    win_ch = channel.of("4") | combine(asm_ch) | create_windows | 
      map { n -> n.get(1) } | splitText() {it.trim() }
    fai_ch = create_windows.out | map { n -> n.get(0) }

    arrow_run_ch = channel.of("4") | combine(asm_ch) | pbmm2_index | combine(pac_ch) | pbmm2_align |
      combine(asm_ch) | combine(fai_ch) | combine(win_ch) | gcc_Arrow 

    if (params.same_specimen) {
      /* calculate peak, and create a genome meryl db */
      peak_ch = channel.of("4")| combine(merylDB_ch) | meryl_peak | splitText() { it.trim() }
      asm_meryl = channel.of(params.k) | combine(asm_ch) | meryl_genome

      /* prepare and run merfin polish */
      newasm_ch = arrow_run_ch | map { n -> [ n.get(0), n.get(2) ] } | groupTuple |
        combineVCF_arrow | reshape_arrow | combine(asm_ch) | combine(asm_meryl) | combine(peak_ch) |
        combine(merylDB_ch) | merfin_polish_arrow | combine(asm_ch) | vcf_to_fasta_arrow
    } else {
      newasm_ch = arrow_run_ch | map { n -> [ n.get(0), n.get(1) ] } | groupTuple |
        merge_consensus
    }
  
  emit:
    newasm_ch
}

process MerquryQV_05 {
    publishDir "${params.outdir}/05_QV/MerquryQV", mode: 'copy'
    publishDir "${params.outdir}/05_QV", mode: 'copy', pattern: "merqury.qv"
    input: tuple path(illumina_db), path(assembly_fasta)
    output: path("*")
    script:
    template 'merquryqv.sh'
}

process bbstat_05 {
    publishDir "${params.outdir}/05_QV/bbstat", mode: 'copy'
    input:  path(assembly_fasta)
    output: path("*")
    script:
    template 'bbstats.sh'
}


// 1st FreeBayes Polish
process align_shortreads {
    publishDir "${params.outdir}/0${i}_FreeBayesPolish/bam", mode: 'symlink'
    input: tuple val(i), path(assembly_fasta), path(illumina_one), path(illumina_two)
    output: tuple val("$i"), path("*.bam"), path("*.bai")
    script:
    template 'align_shortreads.sh'
}
// -m 5G -@ 36

process freebayes {
    publishDir "${params.outdir}/0${i}_FreeBayesPolish/vcf", mode: 'symlink'
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
    publishDir "${params.outdir}/0${i}_FreeBayesPolish/merfin", mode: 'copy'
    input: tuple val(i), path(illumina_reads)
    output: path("*")
    script:
    template 'jellyfish_peak.sh'
}

process meryl_genome_fb {
    publishDir "${params.outdir}/0${i}_FreeBayesPolish/merfin", mode: 'symlink'
    input: tuple val(i), val(k), path(illumina_read)
    output: path("*.meryl")
    script:
    template 'meryl_count.sh'
}

process merfin_polish {
    publishDir "${params.outdir}/0${i}_FreeBayesPolish/merfin", mode: 'symlink'
    input: tuple val(i), path(vcf), path(genome_fasta), path(genome_meryl), val(peak), path(meryldb)
    output: tuple val("$i"), path("*merfin.polish.vcf")
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

workflow FREEBAYES_06 {
  take:
    asm_ch
    ill_ch
    peak_ch
    merylDB_ch
  main:
    win_ch = channel.of("6") | combine(asm_ch) | create_windows | 
      map { n -> n.get(1) } | splitText() {it.trim() }
    fai_ch = create_windows.out | map { n -> n.get(0) }
    asm_meryl = channel.from("6", params.k) | collect | combine(asm_ch) | meryl_genome_fb

    new_asm_ch = channel.of("6") | combine(asm_ch) | combine(ill_ch.collect()) | align_shortreads |
      combine(asm_ch) | combine(fai_ch) | combine(win_ch) |
      freebayes | groupTuple | combineVCF | combine(asm_ch) | 
      combine(asm_meryl) | combine(peak_ch) | combine(merylDB_ch) |  merfin_polish | combine(asm_ch) |
      vcf_to_fasta
  emit:
    new_asm_ch
}

process MerquryQV_07 {
    publishDir "${params.outdir}/07_QV/MerquryQV", mode: 'copy'
    publishDir "${params.outdir}/07_QV", mode: 'copy', pattern: "merqury.qv"
    input: tuple path(illumina_db), path(assembly_fasta)
    output: path("*")
    script:
    template 'merquryqv.sh'
}

process bbstat_07 {
    publishDir "${params.outdir}/07_QV/bbstat", mode: 'copy'
    input:  path(assembly_fasta)
    output: path("*")
    script:
    template 'bbstats.sh'
}

// 2nd Freebayes polish
workflow FREEBAYES_08 {
  take:
    asm_ch
    ill_ch
    peak_ch
    merylDB_ch
  main:
    win_ch = channel.of("8") | combine(asm_ch) | create_windows | 
      map { n -> n.get(1) } | splitText() {it.trim() }
    fai_ch = create_windows.out | map { n -> n.get(0) }
    asm_meryl = channel.from("8", params.k) | collect | combine(asm_ch) | meryl_genome_fb

    new_asm_ch = channel.of("8") | combine(asm_ch) | combine(ill_ch.collect()) | align_shortreads |
      combine(asm_ch) | combine(fai_ch) | combine(win_ch) |
      freebayes | groupTuple | combineVCF | combine(asm_ch) |
      combine(asm_meryl) | combine(peak_ch) | combine(merylDB_ch) |  merfin_polish | combine(asm_ch) |
      vcf_to_fasta
  emit:
    new_asm_ch
}

process MerquryQV_09 {
    publishDir "${params.outdir}/09_QV/MerquryQV", mode: 'copy'
    publishDir "${params.outdir}/09_QV", mode: 'copy', pattern: "merqury.qv"
    input: tuple path(illumina_db), path(assembly_fasta)
    output: path("*")
    script:
    template 'merquryqv.sh'
}

process bbstat_09 {
    publishDir "${params.outdir}/09_QV/bbstat", mode: 'copy'
    input:  path(assembly_fasta)
    output: path("*")
    script:
    template 'bbstats.sh'
}

process addMito {
    publishDir "${params.outdir}/00_Preprocess", mode: 'copy'
    input: tuple path(asm_file), path(mito_file)
    output: path("${asm_file.simpleName}_mito.fasta")
    script:
    """
    #! /usr/bin/env bash
    cat ${asm_file} ${mito_file} > ${asm_file.simpleName}_mito.fasta
    """
}

workflow {
    // Setup input channels, starting assembly (asm), Illumina reads (ill), and pacbio reads (pac)
    asm_ch = channel.fromPath(params.primary_assembly, checkIfExists:true)
    ill_ch = channel.fromFilePairs(params.illumina_reads, checkIfExists:true)
    pac_ch = channel.fromPath(params.pacbio_reads, checkIfExists:true)
    k_ch   = channel.of(params.k) // Either passed in or autodetect (there's a script for this)

    //if(params.mito_assembly){
    //  mito_ch = channel.of(params.mito_assembly, checkIfExists:true)
    //  asm_ch = channel.fromPath(params.primary_assembly, checkIfExists:true) | combine(mito_ch) | addMito
    //} else {
    //  asm_ch = channel.fromPath(params.primary_assembly, checkIfExists:true)
    //}

    // Step 0: Preprocess illumina files from bz2 to gz files
    // Instead of a flag, auto detect, however it must be in the pattern, * will fail
    if(params.illumina_reads =~ /bz2$/){
      pill_ch = ill_ch | bz_to_gz | map { n -> n.get(1) } | flatten
    }else{
      pill_ch = ill_ch | map { n -> n.get(1) } | flatten
    }

    // Step 1: Check quality of assembly with Merqury and length dist. with bbstat
    merylDB_ch = k_ch | combine(pill_ch) | meryl_count | collect | meryl_union 
    merylDB_ch | combine(asm_ch) | MerquryQV_01
    asm_ch | bbstat_01    

    // Step 2: Arrow Polish with PacBio reads
    asm_arrow_ch = ARROW_02(asm_ch, pac_ch)

    // Step 3: Check quality of new assembly with Merqury 
    merylDB_ch | combine(asm_arrow_ch) | MerquryQV_03
      asm_arrow_ch | bbstat_03
 
    // if the primary assembly came from falcon unzip, skip the 2nd arrow polish
    if(!params.falcon_unzip) {
       // Step 4: Arrow Polish with PacBio reads
       asm_arrow2_ch = ARROW_04(asm_arrow_ch, pac_ch, merylDB_ch)
       // Step 5: Check quality of new assembly with Merqury 
       merylDB_ch | combine(asm_arrow2_ch) | MerquryQV_05
       asm_arrow2_ch | bbstat_05
     } else {
       asm_arrow2_ch = asm_arrow_ch
     }
//     asm_arrow2_ch | view
//     
     // Step 6: FreeBayes Polish with Illumina reads
     //peak_ch = channel.of("6") | combine( pill_ch.collect() | map { n-> [n]} ) | jellyfish_peak | splitText() { it.trim()}
     //peak_ch | view
     peak_ch = channel.of("79") // fix this later
     asm_freebayes_ch = FREEBAYES_06(asm_arrow2_ch, pill_ch, peak_ch, merylDB_ch)
     merylDB_ch | combine(asm_freebayes_ch) | MerquryQV_07
     asm_freebayes_ch | bbstat_07
 
     // Step 8: FreeBayes Polish with Illumina reads
     asm_freebayes2_ch = FREEBAYES_08(asm_freebayes_ch, pill_ch, peak_ch, merylDB_ch)
     merylDB_ch | combine(asm_freebayes2_ch) | MerquryQV_09
     asm_freebayes2_ch | bbstat_09
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
