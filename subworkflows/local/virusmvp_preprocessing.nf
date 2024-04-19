#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// import modules
include { EXTRACTVARIANTS       } from '../../modules/local/extractVariants'
include { EXTRACTMETADATA       } from '../../modules/local/extractMetadata'
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

        if(params.variable){
            variable = params.variable
        }
        else{
            variable = []
        }

        criteria = Channel.of(params.grouping_criteria)
        EXTRACTVARIANTS(variants, metadata, variantfile, virusseq, criteria, params.variable, time=true)
        
        EXTRACTVARIANTS.out.txt
            .splitText() 
            .map{id, voc -> tuple([[id:voc.tokenize(':')[1].replaceAll(" ", "_").trim()], voc.replaceAll(" ", "_").trim()])} 
            .set{ ch_group }
        
        
        if(EXTRACTVARIANTS.out.log){
            logfile = EXTRACTVARIANTS.out.log
        }
        EXTRACTMETADATA(metadata, ch_group, time=true)
        
        if (params.grouping_criteria == "time"){
            ids_channel = EXTRACTMETADATA.out.txt.map {it-> it[1]} .flatten()
            ids_channel.view()
            }
        else{
            ids = EXTRACTMETADATA.out.txt
            ids_channel = ids.map{it[1]}
        }
        SEQKIT_GREP(sequences, ids_channel)
        
    emit:
        metadata = EXTRACTMETADATA.out.tsv
        sequences = SEQKIT_GREP.out.filter
        logfile = logfile
}
