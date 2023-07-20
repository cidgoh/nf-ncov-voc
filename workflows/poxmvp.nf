#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// path to required files
/*params.refdb = "$baseDir/assets/virus_referenceGenome"
params.ref_gff = "$baseDir/assets/virus_genomeFeatures"
params.prob_sites = "$baseDir/assets/ncov_problematicSites"
//params.genome_annotation = "$baseDir/assets/ncov_genomeAnnotation/MN908947.3.gff3"
params.functional_annotation = "$baseDir/assets/ncov_functionalAnnotation"
params.gene_coordinates = "$baseDir/assets/virus_refgeneCoordinates"
params.mutation_names = "$baseDir/assets/ncov_multiNames"
params.surveillance_indicators = "$baseDir/assets/ncov_surveillanceIndicators"
*/
// include subworkflows


include {PREPROCESSING          } from '../subworkflows/local/virusmvp_preprocessing'
include {VARIANT_CALLING        } from '../subworkflows/local/virusmvp_variantcalling'
include {ANNOTATION             } from '../subworkflows/local/virusmvp_annotation'
include {surveillance           } from '../subworkflows/local/virusmvp_surveillance'
include {GVF_PROCESSING_ANNOTATION             } from '../subworkflows/local/virusmvp_gvf_processing_annotations'
include {QUALITYCONTROL         } from '../subworkflows/local/virusmvp_qc'
//include { BAM_VARIANT_CALLING_SORT_FREEBAYES_BCFTOOLS } from '../subworkflows/nf-core/bam_variant_calling_sort_freebayes_bcftools/main'


// include modules
include { SNPEFF_BUILD          } from '../modules/local/snpeff_build'
include { PANGOLIN              } from '../modules/nf-core/pangolin/main'
include { EXTRACTVARIANTS       } from '../modules/local/extractVariants'
include { MERGE_PANGOLIN_METADATA } from '../modules/local/merge_pangolin_report'
include { SAMTOOLS_FAIDX       } from '../modules/nf-core/samtools/faidx/main'




workflow POXMVP {

    main:      
        // Building SNPEFF database for COVID
        ch_snpeff_db     = Channel.empty()
        ch_snpeff_config = Channel.empty()
        ch_voc = Channel.empty()
        
        SNPEFF_BUILD (
                params.viral_genome,
                params.viral_gbk
        )
        ch_snpeff_db     = SNPEFF_BUILD.out.db
        ch_snpeff_config = SNPEFF_BUILD.out.config

        SAMTOOLS_FAIDX([[id: params.viral_genome_id], params.viral_genome], [[],[]])
        ch_viral_fai=SAMTOOLS_FAIDX.out.fai


        if(params.meta){
                Channel.fromPath( "$params.meta", checkIfExists: true)
                .set{ ch_metadata }
                }

        if(params.seq){
                Channel.fromPath( "$params.seq", checkIfExists: true)
                .set{ ch_seq }
                }

        seq = file(params.seq, checkIfExists: true)
        id = seq.getSimpleName()
        sequences = [ [ id:id ], [ seq ] ]

        meta = file(params.meta, checkIfExists: true)
        metadata = [ [ id:'Metadata' ], [ meta ] ]

        // Pangolin and extraction of metadata for COVID
        if (!params.skip_pangolin){
            PANGOLIN( sequences )
            MERGE_PANGOLIN_METADATA(metadata, PANGOLIN.out.report)
            metadata = MERGE_PANGOLIN_METADATA.out.tsv
        }
        
        if(params.mode == 'reference'){
            variant = file(params.variant, checkIfExists: true)
            variants = [ [ id:'SARS-CoV-2' ], [ variant ] ]
            EXTRACTVARIANTS(variants, metadata)
            EXTRACTVARIANTS.out.txt
                .splitText()
                .set{ ch_voc }
            }
            

        if(params.mode == 'reference'){      
            PREPROCESSING(metadata, sequences, ch_voc)
            metadata =PREPROCESSING.out.metadata 
            sequences=PREPROCESSING.out.sequences
        }

        sequences
            .map { key, fasta_files ->
              tuple( [[id:fasta_files.getBaseName(2)], [fasta_files]] )
            }
            .set{sequences_grouped}

        
        if(!params.skip_qc){
            sequences
                .map { [it[1]] }
                .collect()
                .map { sequences -> [ [id:"qc"], sequences ] }
                .set { ch_collected_sequences} 

            QUALITYCONTROL(sequences_grouped, ch_collected_sequences)
            sequences_grouped=QUALITYCONTROL.out.sequences_grouped
            ch_stats=QUALITYCONTROL.out.stats
        }
        
        VARIANT_CALLING(sequences_grouped, params.viral_genome, params.viral_genome_fai)
        annotation_vcf=VARIANT_CALLING.out.vcf
        
        ANNOTATION(annotation_vcf, ch_snpeff_db, ch_snpeff_config, params.viral_genome, ch_stats)
        annotation_gvf=ANNOTATION.out.gvf

        GVF_PROCESSING_ANNOTATION(annotation_gvf)


        

        //ch_vcf=variant_calling.out.ch_vcf
        //ch_stats=variant_calling.out.ch_stats
        //ch_metadata=variant_calling.out.ch_metadata
        /*
        annotation(ch_vcf, ch_geneannot, ch_funcannot, ch_genecoord, ch_mutationsplit, ch_variant, ch_stats)
        ch_gvf_surveillance=annotation.out.ch_gvf_surv
        ch_variant=annotation.out.ch_variant
        ch_stats=annotation.out.ch_stats

        surveillance(ch_gvf_surveillance, ch_variant, ch_stats, ch_surveillanceIndicators, ch_metadata )

        }

        
        if (params.mode == 'user' && params.userfile){
                input_file = file(params.userfile)
                ch_metadata=Channel.empty()
                ch_voc=Channel.empty()
            

                if (input_file.getExtension() == "fasta" || input_file.getExtension() == "fa"){
                    Channel.fromPath( "$params.userfile", checkIfExists: true)
                    .set{ ch_seq }

                variant_calling(ch_voc, ch_metadata, ch_seq, ch_ref, ch_refgff, ch_reffai)
                ch_stats=variant_calling.out.ch_stats
                ch_vcf=variant_calling.out.ch_vcf
            
        
                annotation(ch_vcf, ch_geneannot, ch_funcannot, ch_genecoord, ch_mutationsplit, ch_variant, ch_stats)
                ch_gvf_surveillance=annotation.out.ch_gvf_surv
                ch_metadata=ch_mutationsplit
                surveillance(ch_gvf_surveillance, ch_variant , ch_stats, ch_surveillanceIndicators, ch_metadata)
                }

                if (input_file.getExtension() == "vcf" || input_file.getExtension() == "tsv"){
                if(input_file.getExtension() == "vcf"){
                    Channel.fromPath( "$params.userfile", checkIfExists: true)
                    .set{ ch_vcf }
                    }
                
                annotation(ch_vcf, ch_geneannot, ch_funcannot, ch_genecoord, ch_mutationsplit, ch_variant, ch_stats)
                ch_gvf_surveillance=annotation.out.ch_gvf_surv
                surveillance(ch_gvf_surveillance, ch_variant , ch_stats, ch_surveillanceIndicators, ch_metadata)
                    
                }
        }

       

        if (params.mode == 'wastewater' ) {
                //
                // SUBWORKFLOW: Read in samplesheet, validate and stage input files
                //
                
                ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
                wastewater(INPUT_CHECK.out.reads, params.host_genome, params.viral_genome)
                ch_stats=wastewater.out.stats
                annotation(wastewater.out.vcf, ch_stats, ch_snpeff_db, ch_snpeff_config)

        }*/

}