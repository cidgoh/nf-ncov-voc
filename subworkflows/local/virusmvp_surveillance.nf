#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// import modules
include { GVF2TSV } from '../../modules/local/gvf2tsv'
include { TSV2PDF } from '../../modules/local/tsv2pdf'

workflow SURVEILLANCE {
  take:
  ch_gvf
  ch_metadata // This will be a value channel, potentially containing null

  main:
  GVF2TSV(ch_gvf)
  tsv = GVF2TSV.out.surveillancetsv
  surveillance_indicators = Channel
    .fromPath(params.surveillance_indicators)
    .ifEmpty { error("Cannot find viral genome file: ${params.surveillance_indicators}") }
    .map { file -> [[id: "surveillance_indicators"], file] }
    .collect()


  if (params.mode == 'reference') {
    ch_combined = tsv.join(ch_metadata, by: 0)
  }
  else {
    ch_combined = tsv.map { id, gvf -> [id, gvf, []] }
  }
  surveillance_indicators_value = surveillance_indicators.first()
  TSV2PDF(ch_combined, surveillance_indicators_value)

  emit:
  tsv = tsv
  pdf = TSV2PDF.out.surveillance_pdf
}
