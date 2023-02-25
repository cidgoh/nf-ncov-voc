#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// import modules

include { surveillanceRawTsv         } from '../modules/local/custom'
include { surveillancePDF            } from '../modules/local/custom'


workflow surveillance {
    take:
      ch_gvf
      ch_variant
      ch_stats
      ch_surveillanceIndicators
      ch_metadata

    main:
      surveillanceRawTsv(ch_gvf, ch_variant.combine(ch_stats))
      surveillancePDF(surveillanceRawTsv.out.surveillancetsv.flatten(), ch_surveillanceIndicators, ch_metadata)
      ch_surv=surveillancePDF.out.surveillance_pdf
    
    emit:
      ch_surv

}
