#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// import modules

include { NCOVSPLITMUTATIONSGVF                     } from '../../modules/local/splitmutations_gvf'
include { NCOVSPLITMUTATIONSPOKAY             } from '../../modules/local/splitmutations_pokay'
include { FUNCTIONALANNOTATION                  } from '../../modules/local/addFunctionalAnnotation'
include { VARIANTANNOTATION                  } from '../../modules/local/addVariantAnnotation'


workflow GVF_PROCESSING_ANNOTATION {
    take:
        annotation_gvf
        
    main:

        
        if(!params.mpox){
            functional_annotation = file(params.funcannot, checkIfExists: true)
            func = [ [ id:params.viral_genome_id ],  functional_annotation  ]
            
            split_names = file(params.mutationsplit, checkIfExists: true)
            split_tsv = [ [ id:params.viral_genome_id ],  split_names  ]
            
            if(!params.skip_splitting_mutations){
                
                NCOVSPLITMUTATIONSPOKAY(
                    func,
                    split_tsv
                )

                NCOVSPLITMUTATIONSGVF(
                    annotation_gvf,
                    split_tsv
                )
            }
            
            FUNCTIONALANNOTATION(
                NCOVSPLITMUTATIONSGVF.out.gvf,
                NCOVSPLITMUTATIONSPOKAY.out.tsv
            )
            annotation_gvf=FUNCTIONALANNOTATION.out.gvf
            

        }
        variant_annotation = file(params.variant, checkIfExists: true)
        variant_tsv = [ [ id:params.viral_genome_id ],  variant_annotation  ]
    
        VARIANTANNOTATION(
            annotation_gvf,
            variant_tsv,
            true
        )
        annotation_gvf=VARIANTANNOTATION.out.gvf
        

    emit:
        annotation_gvf
        
}
