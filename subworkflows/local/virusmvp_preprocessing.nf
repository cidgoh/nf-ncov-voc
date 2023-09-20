#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// import modules
include { EXTRACTVARIANTS       } from '../../modules/local/extractVariants'
include { extractMetadata       } from '../../modules/local/extractMetadata'
include { SEQKIT_GREP           } from '../../modules/nf-core/seqkit/grep/main'


workflow PREPROCESSING {
    take:
        metadata
        sequences



    main:
        variant = file(params.variant, checkIfExists: true)
        variants = [ [ id:params.viral_genome_id ], variant ]
        if(params.virusseq_update){
                     
            EXTRACTVARIANTS(variants, metadata, [], true)
        }
        else if(!params.skip_variantparsing){
                     
            EXTRACTVARIANTS(variants, metadata, true)
        }
        else{
            EXTRACTVARIANTS(variants, metadata, [])
        }
        
        
        EXTRACTVARIANTS.out.txt
            .splitText() 
            .map{id, voc -> tuple([[id:voc.trim()], voc.trim()])} 
            .set{ ch_voc }

        if(params.grouping_criteria == 'lineage'){
          extractMetadata(metadata, ch_voc, [])
          ids=extractMetadata.out.txt
        }
        else{
          extractMetadata(metadata, [], true)
          //ids=extractMetadata.out.txt.flatten()
          //metadata=extractMetadata.out.tsv.flatten()
        }
        

        SEQKIT_GREP(sequences, ids.map{it[1]})

  emit:
      metadata = extractMetadata.out.tsv
      sequences = SEQKIT_GREP.out.filter
      

}
