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
  surveillance_indicators = Channel.value([id: "surveillance", file: params.surveillance_indicators])
  if (params.mode == 'reference') {
    metadata = ch_metadata
  }
  else {
    metadata = [[], []]
  }
  // Use a conditional operator to handle the case where metadata might be null
  TSV2PDF(tsv, surveillance_indicators, metadata)
}
