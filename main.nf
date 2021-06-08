#! /usr/bin/env nextflow

nextflow.enable.dsl=2

params.primary_assembly="*p_ctg.fasta"
params.illumina_reads="*.fastq.gz"
params.pacbio_reads="*_subreads.fastq.gz"

// 1st Merqury QV value
process meryl_count_01 {
    publishDir "${params.outdir}/01_MerquryQV", mode: 'symlink'
    input: tuple val(k), path(illumina_read)
    output: path("*.meryl")
    shell:
    """
    #! /usr/bin/env bash
    meryl count k=${k} ${illumina_read.simpleName}.meryl ${illumina_read}
    """
}

process meryl_union_01 {
    publishDir "${params.outdir}/01_MerquryQV", mode: 'symlink'
    input: path(illumina_meryls)
    output: path("illumina.meryl")
    shell:
    """
    #! /usr/bin/env bash
    meryl union-sum output illumina.meryl ${illumina_meryls}
    """
}

process MerquryQV_01 {
    publishDir "${params.outdir}/01_MerquryQV", mode: 'symlink'

    input: tuple path(illumina_db), path(assembly_fasta)
    output: path("*")
    shell:
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
    shell:
    """
    #! /usr/bin/env bash
    pbmm2 index ${assembly_fasta} ${assembly_fasta}.mmi
    """
}
process pbmm2_align_01 {
    publishDir "${params.outdir}/02_ArrowPolish", mode: 'symlink'

    input:tuple path(assembly_fasta), path(assembly_mmi), path(pacbio_read)
    output: tuple path("*_aln.bam"), path("*_aln.bai")
    shell:
    """
    #! /usr/bin/env bash
    PROC=`nproc`
    pbmm2 align -j \$PROC ${assembly_fasta} ${pacbio_read} | samtools sort --threads 4 - > ${pacbio_read.simpleName}_aln.bam
    samtools index -@ ${PROC} ${pacbio_read.simpleName}_aln.bam
    """
}

process samtools_faidx_01 {

}

workflow {
    // Setup input channels, starting assembly (asm), Illumina reads (ill), and pacbio reads (pac)
    asm_ch = channel.fromPath(params.primary_assembly, checkIfExists:true)
    ill_ch = channel.fromPath(params.illumina_reads, checkIfExists:true)
    pac_ch = channel.fromPath(params.pacbio_reads, checkIfExists:true)

    // Step 1: Check quality of assembly with Merqury
    channel.of("21") | combine(ill_ch) | meryl_count_01 | collect | meryl_union_01 | combine(asm_ch) | MerquryQV_01

    // Step 2: Arrow Polish
    asm_ch | pbmm2_index_01 | combine(pac_ch) | pbmm2_align_01



}