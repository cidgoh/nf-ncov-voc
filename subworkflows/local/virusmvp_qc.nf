#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// import modules
include { BBMAP        } from '../../modules/local/bbmap_reformat'
include { SEQKIT_STATS } from '../../modules/nf-core/seqkit/stats/main'

workflow QUALITYCONTROL {
    take:
    sequences
    ch_collected_sequences

    main:
    BBMAP(sequences)
    SEQKIT_STATS(ch_collected_sequences)

    emit:
    sequences = BBMAP.out.fasta
    stats     = SEQKIT_STATS.out.stats
}
