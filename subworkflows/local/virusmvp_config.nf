#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// import modules
//include { CONVERTGFFTOJSON          } from '../../modules/local/convertgff2json'
include { SNPEFF_BUILD   } from '../../modules/local/snpeff_build'
include { SAMTOOLS_FAIDX } from '../../modules/nf-core/samtools/faidx/main'

workflow CONFIGURE_VIRUSMVP {
    main:
    ch_snpeff_db = Channel.empty()
    ch_snpeff_config = Channel.empty()
    //ch_viral_fai = Channel.empty()

    SNPEFF_BUILD(
        params.viral_genome,
        params.viral_gbk
    )

    ch_snpeff_db = SNPEFF_BUILD.out.db
    ch_snpeff_config = SNPEFF_BUILD.out.config

    emit:
    ch_snpeff_db
    ch_snpeff_config
}
