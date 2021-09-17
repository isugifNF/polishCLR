#! /usr/bin/env nextflow

nextflow.enable.dsl=2

// === Import Modules
// create meryl database
include { meryl_count; meryl_union; meryl_peak } from './modules/qv.nf'
// QV checks
include { MerquryQV as MerquryQV_01; MerquryQV as MerquryQV_03; MerquryQV as MerquryQV_05; MerquryQV as MerquryQV_07; MerquryQV as MerquryQV_09 } from './modules/qv.nf'
include {bbstat as bbstat_01; bbstat as bbstat_03; bbstat as bbstat_05; bbstat as bbstat_07; bbstat as bbstat_09} from './modules/qv.nf'

// Arrow workflows
include { ARROW as ARROW_02; ARROW_MERFIN as ARROW_04} from './modules/arrow.nf'
// FreeBayes workflows
include { FREEBAYES as FREEBAYES_06; FREEBAYES as FREEBAYES_08 } from './modules/freebayes.nf'
// Preprocess and helper functions
include {bz_to_gz; MERGE_FILE as MERGE_FILE_00; MERGE_FILE_TRIO; SPLIT_FILE as SPLIT_FILE_03; SPLIT_FILE as SPLIT_FILE_09b} from './modules/helper_functions.nf'
// Other
include { PURGE_DUPS as PURGE_DUPS_03b} from './modules/purge_dups.nf'
include { BUSCO } from './modules/busco.nf'


def helpMessage() {
  log.info isuGIFHeader()
  log.info """
   Usage:
   The typical command for running the pipeline are as follows:
   nextflow run main.nf --primary_assembly "*fasta" --illumina_reads "*{1,2}.fastq.bz2" --pacbio_reads "*_subreads.bam" -resume

   Mandatory arguments:
   --illumina_reads               paired end illumina reads, to be used for Merqury QV scores, and freebayes polish primary assembly
   --pacbio_reads                 pacbio reads in bam format, to be used to arrow polish primary assembly
   --mito_assembly                mitocondrial assembly will be concatinated to the assemblies before polishing [default: false]

   Either FALCON (or FALCON Unzip) assembly:
   --primary_assembly             genome assembly fasta file to polish
   --alt_assembly                 if alternate/haplotig assembly file is provided, will be concatinated to the primary assembly before polishing [default: false]
   --falcon_unzip                 if primary assembly has already undergone falcon unzip [default: false]. If true, will Arrow polish once instead of twice.
   
   Or TrioCanu assembly
   --paternal_assembly            paternal genome assembly fasta file to polish
   --maternal_assembly            maternal genome assembly fasta file to polish

   Optional modifiers   
   --species                      if a string is given, rename the final assembly by species name [default:false]
   --k                            kmer to use in MerquryQV scoring [default:21]
   --same_specimen                if illumina and pacbio reads are from the same specimin [default: true].
   --meryldb                      path to a prebuilt meryl database, built from the illumina reads. If not provided, tehen build.
   
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
if ( params.help || !params.illumina_reads || !params.pacbio_reads ) {
  helpMessage()
  exit 0
}

if ( !params.primary_assembly && !params.paternal_assembly ) {
  helpMessage()
  exit 0
}

// If user uses --profile, exit early. The param should be -profile (one hyphen)
if ( params.profile ) {
  helpMessage()
  println("Instead of --profile, use -profile")
  exit 0
}

workflow {
  // === Setup input channels
  // Option 1: read in FALCON assembly
  if( params.alt_assembly ){ // FALCON or FALCON-Unzip
    mito_ch = channel.fromPath(params.mito_assembly, checkIfExists:true)
    alt_ch = channel.fromPath(params.alt_assembly, checkIfExists:true)
    pri_ch = channel.fromPath(params.primary_assembly, checkIfExists:true) |
      map { n -> n.copyTo("${params.outdir}/00_Preprocess/${params.species}_pri.fasta")} 
    
    asm_ch = pri_ch | combine(alt_ch) | combine(mito_ch) | MERGE_FILE_00

  } else if ( params.primary_assembly ) {
    asm_ch = channel.fromPath(params.primary_assembly, checkIfExists:true) |
      map { n -> n.copyTo("${params.outdir}/00_Preprocess/${params.species}_pri.fasta")} 
  }

  // Option 2: read in TrioCanu assembly
  if( params.paternal_assembly ) {
    mito_ch = channel.fromPath(params.mito_assembly, checkIfExists:true)
    pat_ch = channel.fromPath(params.paternal_assembly, checkIfExists:true) | 
      map { n -> n.copyTo("${params.outdir}/00_Preprocess/${params.species}_pat.fasta")}
    mat_ch = channel.fromPath(params.maternal_assembly, checkIfExists:true) |
      map { n -> n.copyTo("${params.outdir}/00_Preprocess/${params.species}_mat.fasta")}

    // should result in Paternal and Material assemblies being polished separately
    asm_ch = pat_ch | concat(mat_ch) | combine(mito_ch) | MERGE_FILE_TRIO
  }
  
  k_ch   = channel.of(params.k) // Either passed in or autodetect (there's a script for this)
  pac_ch = channel.fromPath(params.pacbio_reads, checkIfExists:true)
  // Step 0: Preprocess illumina files from bz2 to gz files. Instead of a flag, auto detect, however it must be in the pattern, * will fail
  if(params.illumina_reads =~ /bz2$/){
    ill_ch = channel.fromFilePairs(params.illumina_reads, checkIfExists:true) | bz_to_gz | map { n -> n.get(1) } | flatten
  }else{
    ill_ch = channel.fromFilePairs(params.illumina_reads, checkIfExists:true) | map { n -> n.get(1) } | flatten
  }
  // Create meryl database and compute peak
  if( params.meryldb ) {  // If prebuilt, will save time
    merylDB_ch = channel.fromPath(params.meryldb, checkIfExists:true)
  } else {
    merylDB_ch = k_ch | combine(ill_ch) | meryl_count | collect | meryl_union
  }
  peak_ch = merylDB_ch | meryl_peak | map { n -> n.get(0) } | splitText() { it.trim() }
  // Step 1: Check quality of assembly with Merqury and length dist. with bbstat   
  merylDB_ch | combine(asm_ch) | MerquryQV_01
  asm_ch | bbstat_01

  if(!params.steptwo) { // TODO: redo this more elegantly later 

    if (!params.falcon_unzip) {
      // Step 2: Arrow Polish with PacBio reads
      asm_arrow_ch = ARROW_02(asm_ch, pac_ch)
      // Step 3: Check quality of new assembly with Merqury 
      merylDB_ch | combine(asm_arrow_ch) | MerquryQV_03
      asm_arrow_ch | bbstat_03
    } else {
      asm_arrow_ch = asm_ch
    }
    
    // purge_dup would go here  
    // (1) split merged into primary, alt, and mito again
    // (2) purge primary, hap merged with alt, purge hap_alt
    // (3) purged primary -> scaffolding pipeline? (might just need a part1, part2 pipeline)
    // (4) merge scaffolded prime, purged alt, and mito
    asm_arrow_ch | SPLIT_FILE_03 |
      map {n -> [n.get(0), n.get(1)] } |
      combine(pac_ch) |
      PURGE_DUPS_03b

    /* BUSCO check will go here */
    PURGE_DUPS_03b.out | map {n -> [n.get(0), n.get(1)] } | flatMap | BUSCO
  } else {
    asm_arrow_ch = asm_ch

    // Step 4: Arrow Polish with PacBio reads
    asm_arrow2_ch = ARROW_04(asm_arrow_ch, pac_ch, peak_ch, merylDB_ch)
    // Step 5: Check quality of new assembly with Merqury 
    merylDB_ch | combine(asm_arrow2_ch) | MerquryQV_05
    asm_arrow2_ch | bbstat_05

    // Step 6: FreeBayes Polish with Illumina reads
    asm_freebayes_ch = FREEBAYES_06(asm_arrow2_ch, ill_ch, peak_ch, merylDB_ch)
    merylDB_ch | combine(asm_freebayes_ch) | MerquryQV_07
    asm_freebayes_ch | bbstat_07
 
    // Step 8: FreeBayes Polish with Illumina reads
    asm_freebayes2_ch = FREEBAYES_08(asm_freebayes_ch, ill_ch, peak_ch, merylDB_ch)
    merylDB_ch | combine(asm_freebayes2_ch) | MerquryQV_09
    asm_freebayes2_ch | bbstat_09

    asm_freebayes2_ch | SPLIT_FILE_09b
    // either incorporate split pat/mat into split_file.sh or create a decision point here
  }
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
