#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// import modules
include {PANGOLIN              } from '../modules/nf-core/modules/pangolin/main'
include {mergePangolinMetadata } from '../modules/local/custom'
include {extractVariants       } from '../modules/local/custom'
include {virusseqMapLineage    } from '../modules/local/custom'


workflow preprocessing {
    take:

      ch_metadata
      ch_seq
      ch_variant


    main:

      if (params.mode == 'reference' && !params.skip_pangolin){
        PANGOLIN( ch_seq )
        mergePangolinMetadata(ch_metadata.combine(PANGOLIN.out.report))
        mergePangolinMetadata.out.lineage_assigned
            .set{ch_metadata}
        extractVariants(ch_variant.combine(mergePangolinMetadata.out.lineage_assigned))
        extractVariants.out.lineages
              .splitText()
              .set{ ch_voc }
        }


      else if (!params.skip_mapping && params.virusseq_meta && params.mode =='reference') {
        Channel.fromPath( "$params.gisaid_metadata", checkIfExists: true)
             .set{ ch_gisaid_metadata }

        virusseqMapLineage(ch_metadata.combine(ch_gisaid_metadata))
        virusseqMapLineage.out.mapped
            .set{ch_metadata}

        extractVariants(ch_variant.combine(virusseqMapLineage.out.mapped))
        extractVariants.out.lineages
            .splitText()
            .set{ ch_voc }
      }

      else if(params.mode =='reference' && params.skip_mapping && params.skip_pangolin) {
        extractVariants(ch_variant.combine(ch_metadata))
        extractVariants.out.lineages
              .splitText()
              .set{ ch_voc }

      }
      else{
        PANGOLIN( ch_seq )
        ch_voc=Channel.empty()
      }

      emit:
      ch_metadata
      ch_voc
      ch_variant

}
