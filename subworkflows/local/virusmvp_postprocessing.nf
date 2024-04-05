#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// import modules
include { GVF_TO_INDEX_LOG       } from '../../modules/local/gvftoindexandlog'
include { MERGE_INDICES         } from '../../modules/local/mergeIndices'
include { MERGE_LOGFILES        } from '../../modules/local/mergeLogfiles'
//include { GVF_TO_INDEX_LOG       } from '../../modules/local/gvftoindexandlog'
//include { POST_LOG       } from '../../modules/local/post_log'
//include { SEQKIT_GREP           } from '../../modules/nf-core/seqkit/grep/main'


workflow POSTPROCESSING {
    take:
        gvf
        logheader

    main:
        GVF_TO_INDEX_LOG(gvf)
        ch_indexfile=Channel.empty()
        if (!params.mutation_indexfile) {
            ch_indexfile = [[],[]]
        }
        else {
            ch_indexfile = [[id:params.virus_accession_id], file(params.mutation_indexfile, checkIfExists: true)]
        }
        GVF_TO_INDEX_LOG.out.index
                .map { [it[1]] }
                .collect()
                .map { indices -> [ [id:"mergingindices"], indices ] }
                .set { ch_indices}
        MERGE_INDICES(ch_indices, ch_indexfile)
        
        GVF_TO_INDEX_LOG.out.log
                .map { [it[1]] }
                .collect()
                .map { logs -> [ [id:"merginglogfiles"], logs ] }
                .set { ch_logs}

        MERGE_LOGFILES(logheader, ch_logs, ch_indexfile)
        //GVF_TO_INDEX_LOG(ch_collected_gvfs, logfile, ch_indexfile)
        //if(!params.skip_posting){
        //    POST_LOG(UPDATE_INDEX_LOG.out.log,  config=file(params.config, checkIfExists: true))
        //}
        
        
}
