#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// import modules

include { GVF2TSV } from '../../modules/local/gvf2tsv'
//include { surveillancePDF } from '../../modules/local/archive/custom'

workflow SURVEILLANCE {
  take:
  ch_gvf

  main:
  GVF2TSV(ch_gvf)
  //surveillancePDF(surveillanceRawTsv.out.surveillancetsv.flatten(), ch_surveillanceIndicators, ch_metadata)
  //ch_surv = surveillancePDF.out.surveillance_pdf
}
