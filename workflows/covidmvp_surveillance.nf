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
      Channel.fromPath( "$baseDir/assets/logo/*.png", checkIfExists: true)
            .set{ ch_logo }
      surveillancePDF(surveillanceRawTsv.out.surveillancetsv.flatten(), ch_surveillanceIndicators.combine(ch_metadata).combine(ch_logo))


}
