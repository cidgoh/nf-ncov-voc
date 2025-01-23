#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// import modules
include { SNPEFF_BUILD   } from '../../modules/local/snpeff_build'
include { SAMTOOLS_FAIDX } from '../../modules/nf-core/samtools/faidx/main'

workflow CONFIGURE_VIRUSMVP {
    main:
    ch_snpeff_db = Channel.empty()
    ch_snpeff_config = Channel.empty()

    // Create a channel for the viral genome file
    ch_viral_genome = Channel
        .fromPath(params.viral_genome)
        .ifEmpty { error("Cannot find viral genome file: ${params.viral_genome}") }

    // Create a channel for the viral GBK file
    ch_viral_gbk = Channel
        .fromPath(params.viral_gbk)
        .ifEmpty { error("Cannot find viral GBK file: ${params.viral_gbk}") }

    SNPEFF_BUILD(
        ch_viral_genome,
        ch_viral_gbk,
    )

    ch_snpeff_db = SNPEFF_BUILD.out.db
    ch_snpeff_config = SNPEFF_BUILD.out.config

    emit:
    ch_snpeff_db
    ch_snpeff_config
}
