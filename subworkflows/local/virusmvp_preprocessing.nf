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
            virusseq=true
        }
        else{
            virusseq=false
        }

        if(params.skip_variantparsing){
            variantfile = []
        }
        else{
            variantfile = true
        }

        criteria = Channel.of(params.grouping_criteria.tokenize(',')).flatten()
        EXTRACTVARIANTS(variants, metadata, variantfile, virusseq, criteria, variable = [], time=true)
        
        EXTRACTVARIANTS.out.txt
            .splitText() 
            .map{id, voc -> tuple([[id:voc.tokenize(':')[1].replaceAll(" ", "_").trim()], voc.replaceAll(" ", "_").trim()])} 
            .set{ ch_group }
        if(EXTRACTVARIANTS.out.log)
        logfile = EXTRACTVARIANTS.out.log
        extractMetadata(metadata, ch_group, time=true)
        ids=extractMetadata.out.txt
        SEQKIT_GREP(sequences, ids.map{it[1]})

    emit:
        metadata = extractMetadata.out.tsv
        sequences = SEQKIT_GREP.out.filter
        logfile = logfile
      
}
