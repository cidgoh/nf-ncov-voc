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


include {PREPROCESSING              } from '../subworkflows/local/virusmvp_preprocessing'
include {VARIANT_CALLING            } from '../subworkflows/local/virusmvp_variantcalling'
include {ANNOTATION                 } from '../subworkflows/local/virusmvp_annotation'
include {surveillance               } from '../subworkflows/local/virusmvp_surveillance'
include {GVF_PROCESSING_ANNOTATION  } from '../subworkflows/local/virusmvp_gvf_processing_annotations'
include {QUALITYCONTROL             } from '../subworkflows/local/virusmvp_qc'
include {VIRALAI                    } from '../subworkflows/local/viralai_datadownload'
//include { BAM_VARIANT_CALLING_SORT_FREEBAYES_BCFTOOLS } from '../subworkflows/nf-core/bam_variant_calling_sort_freebayes_bcftools/main'


// include modules
include { SNPEFF_BUILD              } from '../modules/local/snpeff_build'
include { PANGOLIN                  } from '../modules/nf-core/pangolin/main'
include { MERGE_PANGOLIN_METADATA   } from '../modules/local/merge_pangolin_report'
include { SAMTOOLS_FAIDX            } from '../modules/nf-core/samtools/faidx/main'
include { NEXTCLADE_DATASETGET      } from '../modules/nf-core/nextclade/datasetget/main'
include { NEXTCLADE_RUN             } from '../modules/nf-core/nextclade/run/main'
include { METADATA_HARMONIZER   } from '../modules/local/harmonize_metadata'



workflow COVIDMVP {

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
        
        if(params.viralai){
            VIRALAI()
            metadata=VIRALAI.out.meta
            seq=VIRALAI.out.seq

            seq
                .map { fasta ->
                tuple( [[id:"viralai_seq"], fasta] )
                }
                .set{sequences}      
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
            }
        
            
            if (!params.skip_nextclade){
                dataset = 'sars-cov-2'
                reference = 'MN908947'
                tag = '2023-08-09T12:00:00Z'
                NEXTCLADE_DATASETGET ( dataset, reference, tag ) 
                NEXTCLADE_RUN ( sequences, NEXTCLADE_DATASETGET.out.dataset )
                //MERGE_PANGOLIN_METADATA(metadata, PANGOLIN.out.report)
                //metadata = MERGE_PANGOLIN_METADATA.out.tsv
            }
            
            if (!params.skip_pangolin){
                PANGOLIN( sequences )
                MERGE_PANGOLIN_METADATA(metadata, PANGOLIN.out.report)
                metadata = MERGE_PANGOLIN_METADATA.out.tsv
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
        
        ANNOTATION(annotation_vcf, ch_snpeff_db, ch_snpeff_config, params.viral_genome, ch_stats)
        annotation_gvf=ANNOTATION.out.gvf

        GVF_PROCESSING_ANNOTATION(annotation_gvf)
        

        //surveillance(ch_gvf_surveillance, ch_variant, ch_stats, ch_surveillanceIndicators, ch_metadata )
    
}