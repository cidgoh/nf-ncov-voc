#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// import modules
include { UPDATE_INDEX_LOG       } from '../../modules/local/update_log_index'
include { POST_LOG       } from '../../modules/local/post_log'
//include { SEQKIT_GREP           } from '../../modules/nf-core/seqkit/grep/main'


workflow POSTPROCESSING {
    take:
        ch_collected_gvfs
        logfile

    main:
        if (!params.indexfile) {
            indexfile = file(params.mutation_index, checkIfExists: true)
            ch_indexfile = [[id:params.virus_accession_id], indexfile]
        }

        UPDATE_INDEX_LOG(ch_collected_gvfs, logfile, ch_indexfile)
        if(!params.skip_posting){
            POST_LOG(UPDATE_INDEX_LOG.out.log,  config=file(params.config, checkIfExists: true))
        }
        
        
}
