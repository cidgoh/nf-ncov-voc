#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// import modules
//include { CONVERTGFFTOJSON          } from '../../modules/local/convertgff2json'
include { SNPEFF_BUILD              } from '../../modules/local/snpeff_build'
include { SAMTOOLS_FAIDX            } from '../../modules/nf-core/samtools/faidx/main'


workflow CONFIGURE_VIRUSMVP {

    take:
        //gff
        
    main:
        //ch_json = Channel.empty()
        ch_snpeff_db     = Channel.empty()
        ch_snpeff_config = Channel.empty()
        ch_viral_fai     = Channel.empty()

        // if (params.gene_color) {
        //    color = file(params.gene_color, checkIfExists: true)
        // } else {
        //    color = []
        // }
        // if (params.gene_alias) {
        //     alias = file(params.gene_alias, checkIfExists: true)
        // } else {
        //     alias = []
        // }
        
        // CONVERTGFFTOJSON(gff, color, alias)
        // ch_json = CONVERTGFFTOJSON.out.json
        
        SNPEFF_BUILD (
                params.viral_genome,
                params.viral_gbk
        )

        ch_snpeff_db     = SNPEFF_BUILD.out.db
        ch_snpeff_config = SNPEFF_BUILD.out.config
        
        SAMTOOLS_FAIDX([[id: params.virus_accession_id], params.viral_genome], [[],[]])
        ch_viral_fai=SAMTOOLS_FAIDX.out.fai
        
    emit:
        //ch_json
        ch_snpeff_db
        ch_snpeff_config
        ch_viral_fai
}
