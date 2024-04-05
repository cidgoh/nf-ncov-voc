#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// include subworkflows
include {PREPROCESSING              } from '../subworkflows/local/virusmvp_preprocessing'
include {VARIANT_CALLING            } from '../subworkflows/local/virusmvp_variantcalling'
include {ANNOTATION                 } from '../subworkflows/local/virusmvp_annotation'
include {surveillance               } from '../subworkflows/local/virusmvp_surveillance'
include {GVF_PROCESSING_ANNOTATION  } from '../subworkflows/local/virusmvp_gvf_processing_annotations'
include {QUALITYCONTROL             } from '../subworkflows/local/virusmvp_qc'
include {VIRALAI                    } from '../subworkflows/local/viralai/viralai_datadownload'
include {CLASSIFICATION             } from '../subworkflows/local/virusmvp_classification'
include {POSTPROCESSING             } from '../subworkflows/local/virusmvp_postprocessing'

// include modules
include { METADATA_HARMONIZER       } from '../modules/local/harmonize_metadata'


workflow COVIDMVP {

    take:
        ch_json
        ch_snpeff_db
        ch_snpeff_config
        ch_viral_fai
        
    main:      
        ch_voc = Channel.empty()
        if(params.viralai){
            VIRALAI()
            metadata=VIRALAI.out.meta
            sequences=VIRALAI.out.seq
        }
        else{
            seq = file(params.seq, checkIfExists: true)
            id = seq.getSimpleName()
            sequences = [ [ id:id ], seq ]
        }

        if (params.mode == "reference"){
            if (!params.viralai){
                meta = file(params.meta, checkIfExists: true)
                id = meta.getSimpleName()
                metadata = [ [ id:id ],  meta ]
            }
            
            if (!params.skip_harmonize){
                config = file(params.metadata_config, checkIfExists: true)
                dataset = params.metadata_source
                METADATA_HARMONIZER(metadata, config, dataset)
                metadata=METADATA_HARMONIZER.out.gz
            }
            else{
                metadata=metadata
            }
            if (!params.skip_classification){
                CLASSIFICATION(metadata, sequences)
                metadata=CLASSIFICATION.out.metadata
                sequences=CLASSIFICATION.out.sequences
            }
            PREPROCESSING(metadata, sequences)
            
            metadata =PREPROCESSING.out.metadata 
            sequences=PREPROCESSING.out.sequences
        
            sequences
                .map { key, fasta_files ->
                tuple( [[id:fasta_files.getBaseName(2)], [fasta_files]] )
                }
                .set{sequences_grouped}
        
            sequences
                .map { [it[1]] }
                .collect()
                .map { sequences -> [ [id:"seqkit_stat"], sequences ] }
                .set { ch_collected_sequences}
        
            if(!params.skip_qc){
                QUALITYCONTROL(sequences_grouped, ch_collected_sequences)
                sequences_grouped=QUALITYCONTROL.out.sequences_grouped
                ch_stats=QUALITYCONTROL.out.stats
            }
        }
        else{
            if(!params.skip_qc){
                QUALITYCONTROL(sequences, sequences)
                sequences_grouped=QUALITYCONTROL.out.sequences_grouped
                ch_stats=QUALITYCONTROL.out.stats
            }    
        }
        
        VARIANT_CALLING(sequences_grouped, params.viral_genome, params.viral_genome_fai)
        annotation_vcf=VARIANT_CALLING.out.vcf
        ANNOTATION(annotation_vcf, ch_snpeff_db, ch_snpeff_config, params.viral_genome, ch_stats, ch_json)
        annotation_gvf=ANNOTATION.out.gvf

        GVF_PROCESSING_ANNOTATION(annotation_gvf)
        
        if(!params.skip_postprocessing){
            /*GVF_PROCESSING_ANNOTATION.out.annotation_gvf
                .map { [it[1]] }
                .collect()
                .map { gvfs -> [ [id:"postprocessing"], gvfs ] }
                .set { ch_collected_gvfs}*/
            POSTPROCESSING(annotation_gvf,PREPROCESSING.out.logfile) 
        }
        
        //surveillance(ch_gvf_surveillance, ch_variant, ch_stats, ch_surveillanceIndicators, ch_metadata )
        
}