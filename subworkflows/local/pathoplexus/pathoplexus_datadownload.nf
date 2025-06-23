#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// import modules


include { PATHOPLEXUS_METADATA                  } from '../../../modules/local/pathoplexus/pathoplexus_metadata'
include { PATHOPLEXUS_SEQUENCES                } from '../../../modules/local/pathoplexus/pathoplexus_sequences' 
include { XZ_DECOMPRESS                                  } from '../../../modules/nf-core/xz/decompress/main'

workflow PATHOPLEXUS {

    main:
        
        PATHOPLEXUS_METADATA()
        meta = PATHOPLEXUS_METADATA.out.csv

        PATHOPLEXUS_SEQUENCES()
        seq=PATHOPLEXUS_SEQUENCES.out.xz
        
        seq
            .map { fasta ->
            tuple( [[id:"viralai_seq"], fasta] )
            }
            .set{sequences}
        
        XZ_DECOMPRESS(sequences)
        seq=XZ_DECOMPRESS.out.file
        

        meta
            .map { csv ->
            tuple( [[id:"viralai_meta"], csv] )
            }
            .set{metadata} 

        PROCESS_VIRALAI_METADATA(metadata, alias)
        meta = PROCESS_VIRALAI_METADATA.out.gz
        
    emit:
        meta
        seq
}
