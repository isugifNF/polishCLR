#! /usr/bin/env nextflow

nextflow.enable.dsl=2

params.primary_assembly="*p_ctg.fasta"
params.illumina_reads="*.fastq.gz"
params.pacbio_reads="*_subreads.fastq.gz"
params.outdir="PolishCLR_Results"
params.k="21"


// Unzip any bz2 files
process bz_to_gz {
    publishDir "${params.outdir}/00_Preprocess", mode: 'symlink'
    input:tuple val(readname), path(illumina_reads)
    output: tuple val(readname), path("*.gz")
    script:
    """
    #! /usr/bin/env bash
    PROC=\$((`nproc`-4))
    parallel -j 2 "bzcat {1} | pigz -p \$PROC > {1}.fastq.gz" ::: *.fastq.bz2
    """
}


// 1st Merqury QV value
process meryl_count_01 {

    input: tuple val(k), path(illumina_read)
    output: path("*.meryl")
    script:
    """
    #! /usr/bin/env bash
    meryl count k=${k} output ${illumina_read.simpleName}.meryl ${illumina_read}
    """
}


process meryl_union_01 {
    publishDir "${params.outdir}/01_MerquryQV", mode: 'symlink'
    input: path(illumina_meryls)
    output: path("illumina.meryl")
    script:
    """
    #! /usr/bin/env bash
    meryl union-sum output illumina.meryl ${illumina_meryls}
    """
}

process MerquryQV_01 {
    publishDir "${params.outdir}/01_MerquryQV", mode: 'symlink'

    input: tuple path(illumina_db), path(assembly_fasta)
    output: path("*")
    script:
    """
    #! /usr/bin/env bash
    $MERQURY/merqury.sh $illumina_db $assembly_fasta ${assembly_fasta.simpleName}
    """
}


// 2nd Arrow Polish (if after Falcon unzip)
process pbmm2_index_01 {
    publishDir "${params.outdir}/02_ArrowPolish", mode: 'symlink'

    input:path(assembly_fasta)
    output:tuple path("$assembly_fasta"), path("*.mmi")
    script:
    """
    #! /usr/bin/env bash
    pbmm2 index ${assembly_fasta} ${assembly_fasta}.mmi
    """
}

process pbmm2_align_01 {
    publishDir "${params.outdir}/02_ArrowPolish", mode: 'symlink'

    input:tuple path(assembly_fasta), path(assembly_mmi), path(pacbio_read)
    output: tuple path("*.bam"), path("*.bai")
    script:
    """
    #! /usr/bin/env bash
    PROC=\$((`nproc`-10))
    mkdir tmp
    pbmm2 align -j \$PROC ${assembly_fasta} ${pacbio_read} | samtools sort -T tmp -m 8G --threads 8 - > ${pacbio_read.simpleName}_aln.bam
    samtools index -@ \${PROC} ${pacbio_read.simpleName}_aln.bam
    """
}


process create_windows_01 {
    publishDir "${params.outdir}/02_ArrowPolish", mode: 'symlink'

    input: path(assembly_fasta)
    output: tuple path("*.fai"), path("win.txt")
    shell:
    """
    #! /usr/bin/env bash
    samtools faidx ${assembly_fasta}
    cat ${assembly_fasta}.fai | awk '{print \$1 ":0-" \$2}' > win.txt
    """
}


process gcc_Arrow_01 {
    publishDir "${params.outdir}/02_ArrowPolish", mode: 'symlink'

    input: tuple val(window), path(assembly_fasta), path(assembly_fai), path(pacbio_bam), path(pacbio_bai)
    output: tuple path("*.fasta"), path("*.vcf")
    script:
    """
    #! /usr/bin/env bash
    gcpp --algorithm=arrow \
      -x 10 -X 120 -q 0 \
      -j 24 -w $window \
      -r ${assembly_fasta} ${pacbio_bam} \
      -o ${window.replace(':','_')}.vcf,${window.replace(':','_')}.fasta
    """
}

// add nproc
process merge_consensus_01 {
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
process MerquryQV_02 {
    publishDir "${params.outdir}/03_MerquryQV", mode: 'symlink'

    input: tuple path(illumina_db), path(assembly_fasta)
    output: path("*")
    script:
    """
    #! /usr/bin/env bash
    $MERQURY/merqury.sh $illumina_db $assembly_fasta ${assembly_fasta.simpleName}
    """
}


// 1st FreeBayes Polish
// Pick minimap2 or bwa-mem2 for the aligner (switch to bwa-mem2)

process align_shortreads_01 {
    publishDir "${params.outdir}/04_FreeBayesPolish", mode: 'symlink'

    input: tuple path(assembly_fasta), path(illumina_one), path(illumina_two)
    output: tuple path("*.bam"), path("*.bai")
    script:
    """
    #! /usr/bin/env bash
    PROC=\$((`nproc`-4))
    mkdir tmp
    bwa-mem2 index ${assembly_fasta}
    bwa-mem2 mem -SP -t \${PROC} ${assembly_fasta} ${illumina_one} ${illumina_two} |
      samtools sort -T tmp -m 8G --threads 4 - > ${illumina_one.simpleName}_aln.bam
    samtools index -@ \${PROC} ${illumina_one.simpleName}_aln.bam
    """
}
// -m 5G -@ 36

process create_windows_02 {
    publishDir "${params.outdir}/04_FreeBayesPolish", mode: 'symlink'

    input: path(assembly_fasta)
    output: tuple path("*.fai"), path("win.txt")
    shell:
    """
    #! /usr/bin/env bash
    samtools faidx ${assembly_fasta}
    cat ${assembly_fasta}.fai | awk '{print \$1 ":0-" \$2}' > win.txt
    """
}


process freebayes_01 {
    publishDir "${params.outdir}/04_FreeBayesPolish", mode: 'symlink'
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
      --vcf ${illumina_bam.simpleName}_${window.replace(':','_').replace('|','_')}.vcf \
      --fasta-reference ${assembly_fasta}
    """
}

process combineVCF_01 {
    publishDir "${params.outdir}/04_FreeBayesPolish", mode: 'symlink'
    input: path(vcfs)
    output: path("consensus.vcf")
    script:
    """
    #! /usr/bin/env bash
    cat ${vcfs.get(0)} | grep "^#" > consensus.vcf
    cat ${vcfs} | grep -v "^#" >> consensus.vcf
    """
}

/*
process vcf_to_fasta_01 {
    /// ? where is this step?
}

// 3rd Merqury QV value
process MerquryQV_03 {
    publishDir "${params.outdir}/05_MerquryQV", mode: 'symlink'

    input: tuple path(illumina_db), path(assembly_fasta)
    output: path("*")
    script:
    """
    #! /usr/bin/env bash
    $MERQURY/merqury.sh $illumina_db $assembly_fasta ${assembly_fasta.simpleName}
    """
}

// 2nd FreeBayes
// Pick minimap2 or bwa-mem2 for the aligner
process align_shortreads_02 {
    publishDir "${params.outdir}/06_FreeBayesPolish", mode: 'symlink'

    input: tuple path(assembly_fasta), path(illumina_reads)
    output: path("*.bam")
    script:
    """
    #! /usr/bin/env bash
    PROC=`nproc`
    minimap2 -ax sr -t ${PROC} $assembly_fasta $illumina_reads |
      samtools view -uhS - |
      samtools sort --threads 4 - > ${illumina_reads.get(0).simpleName}.bam
    """
}
// samtools sort needs a tmp directory... add this in later
// -m 5G -@ 36

process freebayes_02 {
    publishDir "${params.outdir}/06_FreeBayesPolish", mode: 'symlink'
    input: tuple val(window), path(assembly_fasta), path(assembly_fai), path(illumina_bam)
    output: tuple path("*.vcf")
    script:
    """
    #! /usr/bin/env bash
    freebayes \
      --region ${window} ${params.options} \
      --bam ${illumina_bam} \
      --vcf ${illumina.simpleName}"_"${window.replace(':','_')}".vcf" \
      --fasta-reference ${assembly_fasta}
    """
}

process combineVCF_02 {
	publishDir "${params.outdir}/06_FreeBayesPolish", mode: 'symlink'
	
	input: path(vcfs)
    output: path("consensus.vcf")

	script:
	"""
    #! /usr/bin/env bash

	cat $vcfs |  vcffirstheader > consensus.vcf 
	"""
}

process vcf_to_fasta_02 {
    /// ? where is this step?
}

// 4th Merqury QV value
process MerquryQV_04 {
    publishDir "${params.outdir}/07_MerquryQV", mode: 'symlink'

    input: tuple path(illumina_db), path(assembly_fasta)
    output: path("*")
    script:
    """
    #! /usr/bin/env bash
    $MERQURY/merqury.sh $illumina_db $assembly_fasta ${assembly_fasta.simpleName}
    """
}
*/

workflow {
    // Setup input channels, starting assembly (asm), Illumina reads (ill), and pacbio reads (pac)
    asm_ch = channel.fromPath(params.primary_assembly, checkIfExists:true)
    ill_ch = channel.fromFilePairs(params.illumina_reads, checkIfExists:true) 
    pac_ch = channel.fromPath(params.pacbio_reads, checkIfExists:true)

    pill_ch = ill_ch | bz_to_gz |map {n -> n.get(1)} | flatten

    // Step 1: Check quality of assembly with Merqury
    channel.of(params.k) | combine(pill_ch) | meryl_count_01 | collect | meryl_union_01 | combine(asm_ch) | MerquryQV_01

 
    // Step 2: Arrow Polish with PacBio reads
    asm_ch | pbmm2_index_01 | combine(pac_ch) | pbmm2_align_01

    fai_ch = asm_ch | create_windows_01 | map { n -> n.get(0) }

    create_windows_01.out | 
      map { n -> n.get(1) } | 
      splitText() {it.trim()} |
      combine(asm_ch) | combine(fai_ch) | combine(pbmm2_align_01.out) | 
      gcc_Arrow_01 |
      map { n -> n.get(0)} |
      collect |
      merge_consensus_01

    asm2_ch = merge_consensus_01.out   // <= New Assembly

    // Step 3: Check quality of new assembly with Merqury (turns out we can reuse the illumina database)
    meryl_union_01.out | combine(asm2_ch) | MerquryQV_02

    // Step 4: FreeBayes Polish with Illumina reads
    asm2_ch | combine(pill_ch.collect()) | align_shortreads_01 | view
    fai2_ch = asm2_ch | create_windows_02 | map { n -> n.get(0) } | view

    // since windows will be the same? (actually double check this...)
    create_windows_02.out | 
      map { n -> n.get(1) } | 
      splitText() {it.trim()} |
      combine(asm2_ch) | combine(fai2_ch) | combine(align_shortreads_01.out) |
      freebayes_01 |
      collect |
      combineVCF_01 // | 
/*      vcf_to_fasta


    asm3_ch = vcf_to_fasta.out         // <= New assembly
  
    // Step 5: Check quality of assembly with Merqury
    meryl_union_01.out | combine(asm3_ch) | MerquryQV_03

    // Step 6: FreeBayes Polish 2nd time
    asm3_ch | combine(ill_ch) | align_shortreads_02
    // since windows will be the same?
    create_windows_01.out | 
      map { n -> n.get(1) } | 
      splitText() {it.trim()} |
      combine(asm3_ch) | combine(fai_ch) | combine(align_shortreads_02.out) |
      freebayes_02 |
      collect |
      combineVCF_02 // | vcf_to_fasta

    asm4_ch = vcf_to_fasta.out         // <= New assembly

    // Step 7: Check quality of assembly with Merqury
    meryl_union_01.out | combine(asm4_ch) | MerquryQV_04
*/    
}