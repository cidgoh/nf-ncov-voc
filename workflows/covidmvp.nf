#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// include subworkflows
include { PREPROCESSING             } from '../subworkflows/local/virusmvp_preprocessing'
include { VARIANT_CALLING           } from '../subworkflows/local/virusmvp_variantcalling'
include { ANNOTATION                } from '../subworkflows/local/virusmvp_annotation'
include { SURVEILLANCE              } from '../subworkflows/local/virusmvp_surveillance'
include { GVF_PROCESSING_ANNOTATION } from '../subworkflows/local/virusmvp_gvf_processing_annotations'
include { QUALITYCONTROL            } from '../subworkflows/local/virusmvp_qc'
include { VIRALAI                   } from '../subworkflows/local/viralai/viralai_datadownload'
include { CLASSIFICATION            } from '../subworkflows/local/virusmvp_classification'
include { POSTPROCESSING            } from '../subworkflows/local/virusmvp_postprocessing'

// include modules
include { METADATA_HARMONIZER       } from '../modules/local/harmonize_metadata'


workflow COVIDMVP {
    take:
    ch_json
    ch_snpeff_db
    ch_snpeff_config

    main:
    //ch_voc = Channel.empty()
    skip_variant_calling = false
    if (params.viralai) {
        VIRALAI()
        metadata = VIRALAI.out.meta
        sequences = VIRALAI.out.seq
    }
    else {
        seq = file(params.seq, checkIfExists: true)
        id = seq.getSimpleName()
        sequences = [[id: id], seq]
    }

    if (params.mode == "reference") {
        if (!params.viralai) {
            meta = file(params.meta, checkIfExists: true)
            id = meta.getSimpleName()
            metadata = [[id: id], meta]
        }

        if (!params.skip_harmonize) {
            config = file(params.metadata_config, checkIfExists: true)
            dataset = params.metadata_source
            METADATA_HARMONIZER(metadata, config, dataset)
            metadata = METADATA_HARMONIZER.out.gz
        }
        else {
            metadata = metadata
        }
        if (!params.skip_classification) {
            CLASSIFICATION(metadata, sequences)
            metadata = CLASSIFICATION.out.merged_metadata
        }

        PREPROCESSING(metadata, sequences)
        metadata = PREPROCESSING.out.metadata
        processed_sequences = PREPROCESSING.out.sequences

        // Step 2: Process SEQKIT_GREP output

        ch_collected_sequences = processed_sequences
            .map { it[1] }
            .collect()
            .map { seqs -> [[id: "seqkit_stat"], seqs] }

        if (!params.skip_qc) {
            QUALITYCONTROL(processed_sequences, ch_collected_sequences)
            sequences_grouped = QUALITYCONTROL.out.sequences
            ch_stats = QUALITYCONTROL.out.stats
        }
    }
    else {
        if (seq.getExtension() == "vcf") {
            skip_variant_calling = true
            annotation_vcf = [[id: seq.getBaseName()], seq]
            ch_stats = []
        }
        else {
            if (!params.skip_qc) {
                QUALITYCONTROL(sequences, sequences)
                sequences_grouped = QUALITYCONTROL.out.sequences
                ch_stats = QUALITYCONTROL.out.stats
            }
        }
    }

    if (!skip_variant_calling == true) {
        VARIANT_CALLING(sequences_grouped, params.viral_genome, params.viral_genome_fai)
        annotation_vcf = VARIANT_CALLING.out.vcf
    }

    ANNOTATION(annotation_vcf, ch_snpeff_db, ch_snpeff_config, params.viral_genome, ch_stats, ch_json)
    annotation_gvf = ANNOTATION.out.gvf

    GVF_PROCESSING_ANNOTATION(annotation_gvf)
    annotated_gvf = GVF_PROCESSING_ANNOTATION.out.annotation_gvf
    if (!params.skip_postprocessing) {
        POSTPROCESSING(annotated_gvf, PREPROCESSING.out.logfile)
    }
    if (!params.skip_surveillance) {
        ch_metadata = params.metadata ? Channel.value(file(params.metadata)) : Channel.value(null)
        SURVEILLANCE(annotated_gvf, ch_metadata)
    }
}
