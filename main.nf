#! /usr/bin/env nextflow

nextflow.enable.dsl=2

// === Import Modules
// Creates Meryl Database and checks QV
include { meryl_count;
          meryl_union;
          meryl_peak;
          MerquryQV as MerquryQV_00;
          MerquryQV as MerquryQV_01;
          MerquryQV as MerquryQV_04;
          MerquryQV as MerquryQV_05;
          MerquryQV as MerquryQV_06;
          bbstat as bbstat_00;
          bbstat as bbstat_01;
          bbstat as bbstat_04;
          bbstat as bbstat_05;
          bbstat as bbstat_06 } from './modules/qv.nf'

// Arrow polishing
include { ARROW as ARROW_02;
          ARROW as ARROW_02b;
          ARROW_MERFIN as ARROW_04;
          ARROW_MERFIN as ARROW_04b } from './modules/arrow.nf'

// FreeBayes polishing
include { FREEBAYES as FREEBAYES_05;
          FREEBAYES as FREEBAYES_05b;
          FREEBAYES as FREEBAYES_06;
          FREEBAYES as FREEBAYES_06b } from './modules/freebayes.nf'

// Preprocess, postprocess, and helper functions
include { bz_to_gz;
          PREFIX_FASTA as PREFIX_PRI;
          PREFIX_FASTA as PREFIX_ALT;
          PREFIX_FASTA as PREFIX_MIT;
          CONCATINATE_FASTA as MERGE_FASTAS_01;
          SPLIT_FASTA as SPLIT_FASTA_01;
          SPLIT_FASTA as SPLIT_FASTA_07;
          bam_to_fasta } from './modules/helper_functions.nf'

// Other
include { PURGE_DUPS as PURGE_DUPS_02;
          PURGE_DUPS_CASE1 as PURGE_DUPS_CASE1 } from './modules/purge_dups.nf'

include { BUSCO } from './modules/busco.nf'

def helpMessage() {
  log.info isuGIFHeader()
  log.info """
   Usage:
   The typical command for running the pipeline are as follows:
   nextflow run main.nf --primary_assembly "*fasta" --illumina_reads "*{1,2}.fastq.bz2" --pacbio_reads "*_subreads.bam" -resume

   Mandatory arguments:
   --illumina_reads               Paired end Illumina reads, to be used for Merqury QV scores, and FreeBayes polish primary assembly
   --pacbio_reads                 PacBio reads in bam format, to be used to Arrow polish primary assembly
   --mitochondrial_assembly       Mitochondrial assembly will be concatenated to the assemblies before polishing [default: false]

   Either FALCON (or FALCON-Unzip) assembly:
   --primary_assembly             Genome assembly fasta file to polish
   --alternate_assembly           If alternate/haplotig assembly file is provided, will be concatenated to the primary assembly before polishing [default: false]
   --falcon_unzip                 If primary assembly has already undergone FALCON-Unzip [default: false]. If true, will Arrow polish once instead of twice.

   Pick Step 1 (Arrow, purgedups) or Step 2 (Arrow, FreeBayes, FreeBayes)
   --step                         Run step 1 or step 2 (default: ${params.step})

   Optional modifiers
   --species                      If a string is given, rename the final assembly by species name [default:false]
   --k                            kmer to use in MerquryQV scoring [default:${params.k}]
   --same_specimen                If Illumina and PacBio reads are from the same specimen [default: true].
   --meryldb                      Path to a prebuilt Meryl database, built from the Illumina reads. If not provided, then build.

   Optional configuration arguments
   --parallel_app                 Link to parallel executable [default: '$parallel_app']
   --bzcat_app                    Link to bzcat executable [default: '$bzcat_app']
   --pigz_app                     Link to pigz executable [default: '$pigz_app']
   --meryl_app                    Link to meryl executable [default: '$meryl_app']
   --merqury_sh                   Link to merqury script [default: '\$MERQURY/merqury.sh']
   --pbmm2_app                    Link to pbmm2 executable [default: '$pbmm2_app']
   --samtools_app                 Link to samtools executable [default: '$samtools_app']
   --gcpp_app                     Link to gcpp executable [default: '$gcpp_app']
   --bwamem2_app                  Link to bwamem2 executable [default: '$bwamem2_app']
   --freebayes_app                Link to freebayes executable [default: '$freebayes_app']
   --bcftools_app                 Link to bcftools executable [default: '$bcftools_app']
   --merfin_app                   Link to merfin executable [default: '$merfin_app']

   Optional parameter arguments
   --parallel_params              Parameters passed to parallel executable [default: '$parallel_params']
   --pbmm2_params                 Parameters passed to pbmm2 align [default: '']
   --minimap2_params              Parameters passed to minimap2 -xmap-pb or -xasm5 [default: '$pbmm2_params']
   --gcpp_params                  Parameters passed to gcpp [default: '$gcpp_params']
   --bwamem2_params               Parameters passed to bwamem2 [default: '$bwamem2_params']
   --freebayes_params             Parameters passed to freebayes [default: '$freebayes_params']
   --purge_dups_params            Parameters passed to purge_dups [default: '$purge_dups_params']
   --busco_params                 Parameters passed to busco [default: '$busco_params']
   --merfin_params                Parameters passed to merfin executable [default: '$merfin_params']


   Optional arguments:
   --outdir                       Output directory to place final output [default: 'PolishCLR_Results']
   --clusterOptions               Cluster options for slurm or sge profiles [default slurm: '-N 1 -n 40 -t 04:00:00'; default sge: ' ']
   --threads                      Number of CPUs to use during each job [default: ${params.threads} ]
   --queueSize                    Maximum number of jobs to be queued [default: ${params.queueSize} ]
   --account                      Some HPCs require you supply an account name for tracking usage.  You can supply that here.
   --help                         This usage statement.
   --check_software               Check if software dependencies are available.

  """
}

// Show help message
if ( params.help) {
  helpMessage()
  exit 0
}

if ( !params.pacbio_reads && !params.check_software ) {
  println("[Missing File(s) Error] polishCLR requires '--pacbio_reads [PacBio bam file]' for ArrowPolish steps.")
  exit 0
}
if ( !params.illumina_reads && !params.check_software ) {
  println("[Missing File(s) Error] polishCLR requires '--illumina_reads [illumina_paired_end_{R1,R2}.fq]' for FreeBayes steps.")
  exit 0
}

if ( !params.primary_assembly && !params.check_software ) {
  println("[Missing File(s) Error] polishCLR requires a '--primary_assembly [primary fasta file]' to polish.")
  exit 0
}

if ( params.profile ) {
  println("[Wrong Flag Error] Instead of --profile, use -profile")
  exit 0
}

def parameters_valid = ['help','monochrome_logs','outdir',
  'primary_assembly','alternate_assembly','mitochondrial_assembly',
  'illumina_reads','pacbio_reads','k','species','meryldb','falcon_unzip','same_specimen','steptwo','step',
  'queueSize','account','threads','clusterOptions',
  'queue-size', 'cluster-options',
  'parallel_app','bzcat_app','pigz_app','meryl_app','merqury_sh','pbmm2_app','minimap2_app','samtools_app',
  'gcpp_app','bwamem2_app','freebayes_app','bcftools_app','merfin_app','pbcstat_app',
  'calcuts_app','split_fa_app','purge_dups_app','get_seqs_app','gzip_app','busco_app','busco_lineage',
  'parallel_params','pbmm2_params','minimap2_params','gcpp_params','bwamem2_params','freebayes_params','merfin_params',
  'purge_dups_params','busco_params','check_software','profile'] as Set

def parameter_diff = params.keySet() - parameters_valid
if (parameter_diff.size() != 0){
   exit 1, "[Pipeline error] Parameter(s) $parameter_diff is(are) not valid in the pipeline!\n"
}

process check_software {
  output: stdout()
  script:
  """
  #! /usr/bin/env bash
  echo "===== Dependencies check ====="

  [[ -z `which $parallel_app` ]]   && echo "${parallel_app}   .... need to install." && ERR=1 || echo "${parallel_app}   .... good. " `${parallel_app} --version | head -n1`
  [[ -z `which $bzcat_app` ]]      && echo "${bzcat_app}      .... need to install." && ERR=1 || echo "${bzcat_app}      .... good. " `${bzcat_app} --help &> temp ; head -n 1 temp`
  [[ -z `which $pigz_app` ]]       && echo "${pigz_app}       .... need to install." && ERR=1 || echo "${pigz_app}       .... good. " `${pigz_app} --version`
  [[ -z `which $meryl_app ` ]]     && echo "${meryl_app}      .... need to install." && ERR=1 || echo "${meryl_app}      .... good. " `${meryl_app} --version &> temp ; head temp`
  [[ -z `which $pbmm2_app` ]]      && echo "${pbmm2_app}      .... need to install." && ERR=1 || echo "${pbmm2_app}      .... good. " `${pbmm2_app} --version`
  [[ -z `which $minimap2_app` ]]   && echo "${minimap2_app}   .... need to install." && ERR=1 || echo "${minimap2_app}   .... good. " `${minimap2_app} --version`
  [[ -z `which $samtools_app` ]]   && echo "${samtools_app}   .... need to install." && ERR=1 || echo "${samtools_app}   .... good. " `${samtools_app} --version | head -n1`
  [[ -z `which $gcpp_app` ]]       && echo "${gcpp_app}       .... need to install." && ERR=1 || echo "${gcpp_app}       .... good. " `${gcpp_app} --version &> temp ; head temp`
  [[ -z `which $bwamem2_app` ]]    && echo "${bwamem2_app}    .... need to install." && ERR=1 || echo "${bwamem2_app}    .... good. " `${bwamem2_app} version`
  [[ -z `which $freebayes_app` ]]  && echo "${freebayes_app}  .... need to install." && ERR=1 || echo "${freebayes_app}  .... good. " `${freebayes_app} --version`
  [[ -z `which $bcftools_app` ]]   && echo "${bcftools_app}   .... need to install." && ERR=1 || echo "${bcftools_app}   .... good. " `${bcftools_app} --version | head -n1`
  [[ -z `which $merfin_app` ]]     && echo "${merfin_app}     .... need to install." && ERR=1 || echo "${merfin_app}     .... good. " `${merfin_app} --version &> temp ; head temp`
  [[ -z `which $pbcstat_app` ]]    && echo "${pbcstat_app}    .... need to install." && ERR=1 || echo "${pbcstat_app}    .... good. " `${pbcstat_app} -h &> temp ; head -n2 temp | tail -n1`
  [[ -z `which $calcuts_app` ]]    && echo "${calcuts_app}    .... need to install." && ERR=1 || echo "${calcuts_app}    .... good. " `${calcuts_app} -h &> temp ; head -n2 temp | tail -n1`
  [[ -z `which $split_fa_app` ]]   && echo "${split_fa_app}   .... need to install." && ERR=1 || echo "${split_fa_app}   .... good. " `${split_fa_app} -h &> temp ; head -n2 temp | tail -n1`
  [[ -z `which $purge_dups_app` ]] && echo "${purge_dups_app} .... need to install." && ERR=1 || echo "${purge_dups_app} .... good. " `${purge_dups_app} -h &> temp ; grep Version temp`
  [[ -z `which $get_seqs_app` ]]   && echo "${get_seqs_app}   .... need to install." && ERR=1 || echo "${get_seqs_app}   .... good. " `${get_seqs_app} -h &> temp ; head -n2 temp | tail -n1`
  [[ -z `which $gzip_app` ]]       && echo "${gzip_app}       .... need to install." && ERR=1 || echo "${gzip_app}       .... good. " `${gzip_app} --version | head -n1`
  [[ -z `which stats.sh` ]]        && echo "stats.sh from bbstat.... need to install." && ERR=1 || echo "stat.sh for bbstat      .... good. " `stats.sh | grep "Usage"`
  [[ -z `which $busco_app ` ]]     && echo "${busco_app}      .... need to install." && ERR=1 || echo "${busco_app}      .... good. " `${busco_app} --version`

  """
}

workflow {
  if ( params.check_software ) {
    check_software()
    | view
  } else {
    // === Setup input channels
    if(params.primary_assembly) {
      primary_assembly_ch = channel.fromPath(params.primary_assembly, checkIfExists:true) 
        | view {file -> "Primary Assembly: $file "}
        | combine(channel.of("pri"))
        | PREFIX_PRI
    } else {
      exit 1, "[Missing File(s) Error] polishCLR requires a '--primary_assembly [primary fasta file]' to polish.\n"
    }
    if(params.alternate_assembly) {
      alternate_assembly_ch = channel.fromPath(params.alternate_assembly, checkIfExists:true) 
        | view {file -> "Alternate Assembly: $file "}
        | combine(channel.of("alt"))
        | PREFIX_ALT
    } else {
      log.info("LOG: Not using alternate assembly since --alternate_assembly [assembly fasta file] was not given.")
      alternate_assembly_ch = channel.empty()
    }
    if(params.mitochondrial_assembly) { // Historically required
      mitochondrial_assembly_ch = channel.fromPath(params.mitochondrial_assembly, checkIfExists:true)
        | view {file -> "Mitochondrial Assembly: $file "}
        | combine(channel.of("mit"))
        | PREFIX_MIT
    } else {
      log.info("LOG: Not using mitochondrial assembly since --mitochondrial_assembly [mitochondrial fasta file] was not given.")
      mitochondrial_assembly_ch = channel.empty()
    }

    assembly_ch = primary_assembly_ch
      | concat(alternate_assembly_ch)  // Will be an empty channel in case 1, will have values in case 2 and 3
      | concat(mitochondrial_assembly_ch)
      | collect
      | map { n -> [n]}
      | combine(channel.of("${params.species}_assembly.fasta"))
      | MERGE_FASTAS_01

    // PacBio and Illumina Reads for polishing
    if(params.pacbio_reads){
      pacbio_reads_ch = channel.fromPath(params.pacbio_reads, checkIfExists:true)
        | view {file -> "PacBio Reads: $file"}
      
      pacbio_all_ch = pacbio_reads_ch
        | collect
        | map {n -> [n]}
    } else {
      exit 1, "[Missing File(s) Error] polishCLR requires '--pacbio_reads [PacBio bam file]' for ArrowPolish steps.\n"
    }

    if(params.illumina_reads){
      // Step 0: Preprocess Illumina files from bz2 to gz files. Instead of a flag, auto detect, however it must be in the pattern, * will fail
      if ( params.illumina_reads =~ /bz2$/ ) {
        illumina_reads_ch = channel.fromFilePairs(params.illumina_reads, checkIfExists:true)
          | view {file -> "Illumina Reads: $file"}
          | bz_to_gz
          | map { n -> n.get(1) }
          | flatten
      } else {
        illumina_reads_ch = channel.fromFilePairs(params.illumina_reads, checkIfExists:true)
          | view {file -> "Illumina Reads: $file"}
          | map { n -> n.get(1) }
          | flatten
      }
    } else {
      exit 1, "[Missing File(s) Error] polishCLR requires '--illumina_reads [illumina_paired_end_{R1,R2}.fq]' for FreeBayes steps.\n"
    }

    k_ch   = channel.of(params.k)
    // Create meryl database and compute peak
    if ( params.meryldb ) {
      merylDB_ch = channel.fromPath(params.meryldb, checkIfExists:true)
    } else {
      log.info("Creating the meryl database from the illumina reads. Passing in a precomputed meryl database via the '--meryldb [path/to/db]' flag may speed up QV checks in subsequent runs.")
      merylDB_ch = k_ch
        | combine(illumina_reads_ch)
        | meryl_count
        | collect
        | meryl_union
    }

    peak_ch = merylDB_ch
      | meryl_peak
      | map { n -> n.get(0) }
      | splitText() { it.trim() }

    // Step 1: Check quality of assembly with Merqury and length dist. with bbstat
    assembly_ch
      | combine(merylDB_ch)
      | MerquryQV_00
    
    assembly_ch
      | bbstat_00

    if ( params.step == 1 ) {

      if (!params.falcon_unzip) {
        // Step 2: Arrow Polish with PacBio reads
        if ( params.primary_assembly ) {
          asm_arrow_ch = ARROW_02(
              channel.of("Step_1/01_ArrowPolish"),
              assembly_ch,
              pacbio_all_ch
            )
        }
        // Step 3: Check quality of new assembly with Merqury
        asm_arrow_ch
          | combine(merylDB_ch)
          | MerquryQV_01

        asm_arrow_ch
          | bbstat_01
      } else {
        asm_arrow_ch = assembly_ch
      }

      pacbio_fasta_ch = pacbio_reads_ch
        | bam_to_fasta
        | collect
        | map {n -> [n]}

      pri_asm_arrow_ch = asm_arrow_ch
        | SPLIT_FASTA_01
        | flatten
        | filter { it.getName() =~ /^pri/ }

      alt_asm_arrow_ch = SPLIT_FASTA_01.out
        | flatten
        | filter { it.getName() =~ /^alt/ }

	    // Case 1 - primary and mito assembly, no alternate.
	    if ( !params.alternate_assembly ) {
        purged_asm_arrow_ch = pri_asm_arrow_ch
          | combine(pacbio_fasta_ch)
          | PURGE_DUPS_CASE1
          | map {n -> n.get(0) } // Only purged primary
	    } else {
        purged_asm_arrow_ch = pri_asm_arrow_ch
          | combine(alt_asm_arrow_ch)
          | combine(pacbio_fasta_ch)
          | PURGE_DUPS_02
          | map {n -> [n.get(0), n.get(1)] } // purged primary and alternative
          | flatten
      }
      purged_asm_arrow_ch
        | BUSCO
      
    } else {
      asm_arrow_ch = assembly_ch

      // Step 4: Arrow Polish with PacBio reads
      asm_arrow2_ch = ARROW_04(
        channel.of("Step_2/04_ArrowPolish"), 
        asm_arrow_ch, 
        pacbio_all_ch,
        peak_ch, 
        merylDB_ch)

      // Step 5: Check quality of new assembly with Merqury
      asm_arrow2_ch
        | combine(merylDB_ch)
        | MerquryQV_04
      
      asm_arrow2_ch
        | bbstat_04

      // Step 6: FreeBayes polish with Illumina reads
      asm_freebayes_ch = FREEBAYES_05(
        channel.of("Step_2/05_FreeBayesPolish"),
        asm_arrow2_ch, 
        illumina_reads_ch, 
        peak_ch, 
        merylDB_ch)

      asm_freebayes_ch
        | combine(merylDB_ch)
        | MerquryQV_05
      
      asm_freebayes_ch
        | bbstat_05 

      // Step 8: FreeBayes polish with Illumina reads
      asm_freebayes2_ch = FREEBAYES_06(
        channel.of("Step_2/06_FreeBayesPolish"), 
        asm_freebayes_ch, 
        illumina_reads_ch, 
        peak_ch, 
        merylDB_ch)

      asm_freebayes2_ch
        | combine(merylDB_ch)
        | MerquryQV_06
      
      asm_freebayes2_ch
        | bbstat_06 

      asm_freebayes2_ch
        | SPLIT_FASTA_07
    }
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
