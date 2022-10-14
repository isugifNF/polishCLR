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
          MERGE_FILE as MERGE_FILE_00;
          MERGE_FILE_TRIO;
          SPLIT_FILE as SPLIT_FILE_02;
          SPLIT_FILE as SPLIT_FILE_07;
          SPLIT_FILE_CASE1 as SPLIT_FILE_CASE1;
          SPLIT_FILE_p as SPLIT_FILE_02p;
          SPLIT_FILE_m as SPLIT_FILE_02m;
          SPLIT_FILE_p as SPLIT_FILE_07p;
          SPLIT_FILE_m as SPLIT_FILE_07m;
          RENAME_FILE as RENAME_PRIMARY;
          RENAME_FILE as RENAME_PAT;
          RENAME_FILE as RENAME_MAT;
          MERGE_FILE_CASE1 as MERGE_CASE1;
          bam_to_fasta } from './modules/helper_functions.nf'

// Other
include { PURGE_DUPS as PURGE_DUPS_02;
          PURGE_DUPS_CASE1 as PURGE_DUPS_CASE1;
          PURGE_DUPS_TRIO as PURGE_DUPS_TRIOp;
          PURGE_DUPS_TRIO as PURGE_DUPS_TRIOm } from './modules/purge_dups.nf'

include { BUSCO as BUSCO_CASE1;
          BUSCO as BUSCO_mat;
          BUSCO } from './modules/busco.nf'

def helpMessage() {
  log.info isuGIFHeader()
  log.info """
   Usage:
   The typical command for running the pipeline are as follows:
   nextflow run main.nf --primary_assembly "*fasta" --illumina_reads "*{1,2}.fastq.bz2" --pacbio_reads "*_subreads.bam" -resume

   Mandatory arguments:
   --illumina_reads               Paired end Illumina reads, to be used for Merqury QV scores, and FreeBayes polish primary assembly
   --pacbio_reads                 PacBio reads in bam format, to be used to Arrow polish primary assembly
   --mitochondrial_assembly       Mitocondrial assembly will be concatinated to the assemblies before polishing [default: false]

   Either FALCON (or FALCON-Unzip) assembly:
   --primary_assembly             Genome assembly fasta file to polish
   --alternate_assembly           If alternate/haplotig assembly file is provided, will be concatinated to the primary assembly before polishing [default: false]
   --falcon_unzip                 If primary assembly has already undergone FALCON-Unzip [default: false]. If true, will Arrow polish once instead of twice.

   Or TrioCanu assembly
   --paternal_assembly            Paternal genome assembly fasta file to polish
   --maternal_assembly            Maternal genome assembly fasta file to polish

   Pick Step 1 (arrow, purgedups) or Step 2 (Arrow, FreeBayes, FreeBayes)
   --step                         Run step 1 or step 2 (default: 1)

   Optional modifiers
   --species                      If a string is given, rename the final assembly by species name [default:false]
   --k                            kmer to use in MerquryQV scoring [default:21]
   --same_specimen                If Illumina and PacBio reads are from the same specimin [default: true].
   --meryldb                      Path to a prebuilt Meryl database, built from the Illumina reads. If not provided, then build.

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

   Optional parameter arguments
   --parallel_params              Parameters passed to parallel executable [default: ' -j 2 ']
   --pbmm2_params                 Parameters passed to pbmm2 align [default: '']
   --minimap2_params              Parameters passed to minimap2 -xmap-pb or -xasm5 [default: '']
   --gcpp_params                  Parameters passed to gcpp [default: ' -x 10 -X 120 -q 0 ']
   --bwamem2_params               Parameters passed to bwamem2 [default: ' -SP ']
   --freebayes_params             Parameters passed to freebayes [default: ' --min-mapping-quality 0 --min-coverage 3 --min-supporting-allele-qsum 0  --ploidy 2 --min-alternate-fraction 0.2 --max-complex-gap 0 ']
   --purge_dups_params            Parameters passed to purge_dups [default: ' -2 -T p_cufoffs ']
   --busco_params                 Parameters passed to busco [default: ' -l insecta_odb10 -m genome -f ']


   Optional arguments:
   --outdir                       Output directory to place final output [default: 'PolishCLR_Results']
   --clusterOptions               Cluster options for slurm or sge profiles [default slurm: '-N 1 -n 40 -t 04:00:00'; default sge: ' ']
   --threads                      Number of CPUs to use during each job [default: 40]
   --queueSize                    Maximum number of jobs to be queued [default: 50]
   --account                      Some HPCs require you supply an account name for tracking usage.  You can supply that here.
   --help                         This usage statement.
   --check_software               Check if software dependencies are available.
  """
}

// Show help message
if ( ( params.help || !params.illumina_reads || !params.pacbio_reads ) && !params.check_software ) {
  log.info("$params.check_software")
  log.info("$params.help || !$params.illumina_reads || !$params.pacbio_reads")
  helpMessage()
  exit 0
}

if ( (!params.primary_assembly && !params.paternal_assembly) && !params.check_software ) {
  helpMessage()
  exit 0
}

// If user uses --profile, exit early. The param should be -profile (one hyphen)
if ( params.profile ) {
  helpMessage()
  println("Instead of --profile, use -profile")
  exit 0
}

def parameters_valid = ['help','monochrome_logs','outdir',
  'primary_assembly','alternate_assembly','paternal_assembly','maternal_assembly','mitochondrial_assembly',
  'illumina_reads','pacbio_reads','k','species','meryldb','falcon_unzip','same_specimen','steptwo','step',
  'queueSize','account','threads','clusterOptions',
  'queue-size', 'cluster-options',
  'parallel_app','bzcat_app','pigz_app','meryl_app','merqury_sh','pbmm2_app','minimap2_app','samtools_app',
  'gcpp_app','bwamem2_app','freebayes_app','bcftools_app','merfin_app','pbcstat_app','hist_plot_py',
  'calcuts_app','split_fa_app','purge_dups_app','get_seqs_app','gzip_app','busco_app','busco_lineage',
  'parallel_params','pbmm2_params','minimap2_params','gcpp_params','bwamem2_params','freebayes_params',
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
  [[ -z `which $busco_app ` ]]     && echo "${busco_app}      .... need to install." && ERR=1 || echo "${busco_app}      .... good. " `${busco_app} --version`

  """
}

workflow {
  if ( params.check_software ) {
    check_software()
    | view
  } else {
    // === Setup input channels
    // Case 2 or Case 3: Primary, alternate, and mito exist
    if ( params.alternate_assembly ) {
      mitochondrial_ch = channel.fromPath(params.mitochondrial_assembly, checkIfExists:true)
      alternate_assembly_ch = channel.fromPath(params.alternate_assembly, checkIfExists:true)
      primary_assembly_ch = channel.fromPath(params.primary_assembly, checkIfExists:true)
        | combine(channel.of("${params.species}_pri.fasta"))
        | RENAME_PRIMARY

      assembly_ch = primary_assembly_ch
        | combine(alternate_assembly_ch)
        | combine(mitochondrial_ch)
        | MERGE_FILE_00

    } else if ( params.primary_assembly ) { // Option 1b: Canu, same as FALCON but without alternative assembly
      if ( params.mitochondrial_assembly ) {  // Merge mitochondrial assembly if available
        mitochondrial_ch = channel.fromPath(params.mitochondrial_assembly, checkIfExists:true)
        primary_assembly_ch = channel.fromPath(params.primary_assembly, checkIfExists:true)
          | combine(channel.of("${params.species}_pri.fasta"))
          | RENAME_PRIMARY

        assembly_ch = primary_assembly_ch
          | combine(mitochondrial_ch)
          | MERGE_CASE1

      } else { // If mitochondrial assembly not available, may still have issues with split_file
        assembly_ch = channel.fromPath(params.primary_assembly, checkIfExists:true)
          | combine(channel.of("${params.species}_pri.fasta"))
          | RENAME_PRIMARY
      }
    } else if ( params.paternal_assembly ) {  // Option 2: read in TrioCanu assembly
      mitochondrial_ch = channel.fromPath(params.mitochondrial_assembly, checkIfExists:true)
      paternal_assembly_ch = channel.fromPath(params.paternal_assembly, checkIfExists:true)
        | combine(channel.of("${params.species}_pat.fasta"))
        | RENAME_PAT
      maternal_assembly_ch = channel.fromPath(params.maternal_assembly, checkIfExists:true)
        | combine(channel.of("${params.species}_mat.fasta"))
        | RENAME_MAT

      // Should result in Paternal and Material assemblies being polished separately
      assembly_ch = paternal_assembly_ch
        | concat(maternal_assembly_ch)
        | combine(mitochondrial_ch)
        | MERGE_FILE_TRIO
    }

    k_ch   = channel.of(params.k)
    pacbio_reads_ch = channel.fromPath(params.pacbio_reads, checkIfExists:true)
    pacbio_all_ch = pacbio_reads_ch
      | collect
      | map {n -> [n]}

    // Step 0: Preprocess Illumina files from bz2 to gz files. Instead of a flag, auto detect, however it must be in the pattern, * will fail
    if ( params.illumina_reads =~ /bz2$/ ) {
      illumina_reads_ch = channel.fromFilePairs(params.illumina_reads, checkIfExists:true)
        | bz_to_gz
        | map { n -> n.get(1) }
        | flatten
    } else {
      illumina_reads_ch = channel.fromFilePairs(params.illumina_reads, checkIfExists:true)
        | map { n -> n.get(1) }
        | flatten
    }
    // Create meryl database and compute peak
    if ( params.meryldb ) {
      merylDB_ch = channel.fromPath(params.meryldb, checkIfExists:true)
    } else {
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
    channel.of("00_Preprocess/00_QV")
      | combine(merylDB_ch)
      | combine(assembly_ch)
      | MerquryQV_00
    channel.of("00_Preprocess/00_bbstat")
      | combine(assembly_ch)
      | bbstat_00

    if ( params.step == 1 ) {

      if (!params.falcon_unzip) {
        // Step 2: Arrow Polish with PacBio reads
        if ( params.primary_assembly ) {
          asm_arrow_ch = ARROW_02(channel.of("Step_1/01_ArrowPolish"), assembly_ch, pacbio_all_ch)
        } else if ( params.paternal_assembly ) {
          asm_arrow_pat_ch = ARROW_02(channel.of("Step_1/01_ArrowPolish_pat"), assembly_ch.first(), pacbio_all_ch)
          asm_arrow_mat_ch = ARROW_02b(channel.of("Step_1/01_ArrowPolish_mat"), assembly_ch.last(), pacbio_all_ch)
          asm_arrow_ch = asm_arrow_pat_ch
            | concat(asm_arrow_mat_ch)
        }

        // Step 3: Check quality of new assembly with Merqury
        channel.of("Step_1/01_QV")
          | combine(merylDB_ch)
          | combine(asm_arrow_ch)
          | MerquryQV_01
        channel.of("Step_1/01_bbstat")
          | combine(asm_arrow_ch)
          | bbstat_01
      } else {
        asm_arrow_ch = assembly_ch
      }

      pacbio_fasta_ch = pacbio_reads_ch
        | bam_to_fasta
        | collect
        | map {n -> [n]}

      if ( params.primary_assembly ) {
	      // Case 1 - primary and mito assembly, no alternate.
	      if ( !params.alternate_assembly ) {
            tmp_ch = channel.of("Step_1/01_ArrowPolish")
              | combine(asm_arrow_ch)
              | SPLIT_FILE_CASE1
              | map {n -> [n.get(0)] }

            channel.of("Step_1/02_Purge_Dups")
              | combine(tmp_ch)
              | combine(pacbio_fasta_ch)
              | PURGE_DUPS_CASE1

            PURGE_DUPS_CASE1.out
              | map {n -> [n.get(0), n.get(1)] }
              | flatMap
              | combine(channel.of("Step_1/02_BUSCO"))
              | map {n -> [n.get(1), n.get(0)]}
              | BUSCO_CASE1
	      }

	      if ( params.alternate_assembly ) {
	        tmp_ch = channel.of("Step_1/01_ArrowPolish")
            | combine(asm_arrow_ch)
            | SPLIT_FILE_02
            | map {n -> [n.get(0), n.get(1)] }
          channel.of("Step_1/02_Purge_Dups")
            | combine(tmp_ch)
            | combine(pacbio_fasta_ch)
            | PURGE_DUPS_02

          PURGE_DUPS_02.out
            | map {n -> [n.get(0), n.get(1)] }
            | flatMap
            | combine(channel.of("Step_1/02_BUSCO"))
            | map {n -> [n.get(1), n.get(0)]}
            | BUSCO

        }

      } else if ( params.paternal_assembly ) {
        // Paternal version
        tmp_ch = asm_arrow_ch.first()
          | SPLIT_FILE_02p
          | map {n -> n.get(0) }
        channel.of("Step_1/02_Purge_Dups_pat")
          | combine(tmp_ch)
          | combine(pacbio_fasta_ch)
          | PURGE_DUPS_TRIOp

        PURGE_DUPS_TRIOp.out
          | map {n -> [n.get(0)] }
          | flatMap
          | combine(channel.of("Step_1/02_BUSCO_pat"))
          | map {n -> [n.get(1), n.get(0)]}
          | combine(channel.of("Step_1/02_BUSCO_pat"))
          | map {n -> [n.get(1), n.get(0)]}
          | BUSCO

        // Maternal version
        tmpm_ch = asm_arrow_ch.last()
          | SPLIT_FILE_02m
          | map {n -> n.get(0) }
        channel.of("Step_1/02_Purge_Dups_mat")
          | combine(tmpm_ch)
          | combine(pacbio_fasta_ch)
          | PURGE_DUPS_TRIOm

        PURGE_DUPS_TRIOm.out
          | map {n -> [n.get(0)] }
          | flatMap
          | combine(channel.of("Step_1/02_BUSCO_mat"))
          | map {n -> [n.get(1), n.get(0)]}
          | BUSCO_mat
      }
    } else {
      asm_arrow_ch = assembly_ch

      // Step 4: Arrow Polish with PacBio reads
      if ( params.primary_assembly ) {
        asm_arrow2_ch = ARROW_04(channel.of("Step_2/04_ArrowPolish"), asm_arrow_ch, pacbio_all_ch, peak_ch, merylDB_ch)
      } else if ( params.paternal_assembly ) {
        asm_arrow2_pat_ch = ARROW_04(channel.of("Step_2/04_ArrowPolish_pat"), asm_arrow_ch.first() , pacbio_all_ch, peak_ch, merylDB_ch)
        asm_arrow2_mat_ch = ARROW_04b(channel.of("Step_2/04_ArrowPolish_mat"), asm_arrow_ch.last(), pacbio_all_ch, peak_ch, merylDB_ch)
        asm_arrow2_ch = asm_arrow2_pat_ch
          | concat(asm_arrow2_mat_ch)
      }
      // Step 5: Check quality of new assembly with Merqury
      channel.of("Step_2/04_QV")
        | combine(merylDB_ch)
        | combine(asm_arrow2_ch)
        | MerquryQV_04
      channel.of("Step_2/04_bbstat")
        | combine(asm_arrow2_ch)
        | bbstat_04

      // Step 6: FreeBayes polish with Illumina reads
      if ( params.primary_assembly ) {
        asm_freebayes_ch = FREEBAYES_05(channel.of("Step_2/05_FreeBayesPolish"), asm_arrow2_ch, illumina_reads_ch, peak_ch, merylDB_ch)
      } else if ( params.paternal_assembly ) {
        asm_freebayes_pat_ch = FREEBAYES_05(channel.of("Step_2/05_FreeBayesPolish_pat"), asm_arrow2_ch.first(), illumina_reads_ch, peak_ch, merylDB_ch)
        asm_freebayes_mat_ch = FREEBAYES_05b(channel.of("Step_2/05_FreeBayesPolish_mat"), asm_arrow2_ch.last(), illumina_reads_ch, peak_ch, merylDB_ch)
        asm_freebayes_ch = asm_freebayes_pat_ch
          | concat(asm_freebayes_mat_ch )
      }

      channel.of("Step_2/05_QV")
        | combine(merylDB_ch)
        | combine(asm_freebayes_ch)
        | MerquryQV_05
      channel.of("Step_2/05_bbstat")
        | combine(asm_freebayes_ch)
        | bbstat_05

      // Step 8: FreeBayes polish with Illumina reads
      if ( params.primary_assembly ) {
        asm_freebayes2_ch = FREEBAYES_06(channel.of("Step_2/06_FreeBayesPolish"), asm_freebayes_ch, illumina_reads_ch, peak_ch, merylDB_ch)
      } else if ( params.paternal_assembly ) {
        asm_freebayes2_pat_ch = FREEBAYES_06(channel.of("Step_2/06_FreeBayesPolish_pat"), asm_freebayes_ch.first(), illumina_reads_ch, peak_ch, merylDB_ch)
        asm_freebayes2_mat_ch = FREEBAYES_06b(channel.of("Step_2/06_FreeBayesPolish_mat"), asm_freebayes_ch.last(), illumina_reads_ch, peak_ch, merylDB_ch)
        asm_freebayes2_ch = asm_freebayes2_pat_ch
          | concat(asm_freebayes2_mat_ch)
      }

      channel.of("Step_2/06_QV")
        | combine(merylDB_ch)
        | combine(asm_freebayes2_ch)
        | MerquryQV_06
      channel.of("Step_2/06_bbstat")
        | combine(asm_freebayes2_ch)
        | bbstat_06

      if (params.primary_assembly) {
        channel.of("Step_2/06_FreeBayesPolish")
          | combine(asm_freebayes2_ch)
          | SPLIT_FILE_07
      } else if ( params.paternal_assembly ) {
        asm_freebayes2_ch.first()
          | SPLIT_FILE_07p
        asm_freebayes2_ch.last()
          | SPLIT_FILE_07m
      }
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
