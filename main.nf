#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// path to required files
params.refdb = "$baseDir/assets/ncov_refdb"
params.ref_gff = "$baseDir/assets/ncov_genomeFeatures"
params.prob_sites = "$baseDir/assets/ncov_problematicSites"
//params.genome_annotation = "$baseDir/assets/ncov_genomeAnnotation/MN908947.3.gff3"
params.functional_annotation = "$baseDir/assets/ncov_functionalAnnotation"
params.gene_coordinates = "$baseDir/assets/ncov_geneCoordinates"
params.mutation_names = "$baseDir/assets/ncov_multiNames"
params.surveillance_indicators = "$baseDir/assets/ncov_surveillanceIndicators"

script_files = "$baseDir/bin"

// include modules
include {printHelp              } from './modules/local/help'
include {cidgohHeader           } from './modules/local/header'
include {workflowHeader         } from './modules/local/wf_header'
include {SNPEFF_BUILD           } from './modules/local/snpeff_build'




// import workflows
include {INPUT_CHECK            } from './subworkflows/local/input_check'
include {preprocessing          } from './workflows/covidmvp_preprocessing'
include {variant_calling        } from './workflows/covidmvp_variantcalling'
include {annotation             } from './workflows/covidmvp_annotation'
include {surveillance           } from './workflows/covidmvp_surveillance'
include {wastewater             } from './workflows/covidmvp_wastewater'




if (params.help){
    log.info cidgohHeader()
    log.info workflowHeader()
    printHelp()
    exit 0
}

if (params.profile){
    println("Profile should have a single dash: -profile")
    System.exit(1)
}

if ( ! params.prefix ) {
     println("Please supply a prefix for your output files with --prefix")
     println("Use --help to print help")
     System.exit(1)
} else {
     if ( params.prefix =~ /\// ){
         println("The --prefix that you supplied contains a \"/\", please replace it with another character")
         System.exit(1)
     }
}

if ( params.mode == 'user' && ! params.userfile ) {
    println("When --mode user, userfile (.vcf or .fasta or .tsv) should e provided with --userfile")
    System.exit(1)
}


// main workflow
workflow {

   main:
      log.info cidgohHeader()
      log.info workflowHeader()

      if (!params.skip_permissions) {
            allFiles = listOfFiles = file("$script_files/*")
            for( def file : allFiles ) {
                  file.setPermissions('rwxr-x--x')
            }
      }
      
      
      ch_snpeff_db     = Channel.empty()
      ch_snpeff_config = Channel.empty()
      SNPEFF_BUILD (
            params.viral_genome,
            params.viral_gbk
      )
      ch_snpeff_db     = SNPEFF_BUILD.out.db
      ch_snpeff_config = SNPEFF_BUILD.out.config

      // Channel.fromPath( "$params.refdb/*.fai", checkIfExists: true)
      //       .set{ ch_reffai }

      // Channel.fromPath( "$params.refdb/MN*.fasta", checkIfExists: true)
      //       .set{ ch_ref }

      // Channel.fromPath( "$params.functional_annotation/*.tsv", checkIfExists: true)
      //       .set{ ch_funcannot }

      // Channel.fromPath( "$params.gene_coordinates/*.json", checkIfExists: true)
      //       .set{ ch_genecoord }

      // Channel.fromPath( "$params.mutation_names/*.tsv", checkIfExists: true)
      //       .set{ ch_mutationsplit }

      // Channel.fromPath( "$params.surveillance_indicators/*.tsv", checkIfExists: true)
      //       .set{ ch_surveillanceIndicators}

      // if(params.variants){
      //       Channel.fromPath( "$params.variants", checkIfExists: true)
      //       .set{ ch_variants }
      // }
      // else{
      //       Channel.fromPath( "$baseDir/assets/ncov_variants/*.tsv", checkIfExists: true)
      //       .set{ ch_variant }
      // }
      
      
      

      if (params.mode == 'wastewater' ) {
            //
            // SUBWORKFLOW: Read in samplesheet, validate and stage input files
            //
            if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
            ch_versions = Channel.empty()
            INPUT_CHECK (
                  ch_input
            )
            ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
            wastewater(INPUT_CHECK.out.reads, params.host_genome, params.viral_genome)
            ch_stats=wastewater.out.stats
            annotation(wastewater.out.vcf, ch_stats, ch_snpeff_db, ch_snpeff_config)

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

      if(params.mode == 'reference'){

        
            if(params.seq){
                  Channel.fromPath( "$params.seq", checkIfExists: true)
                  .set{ ch_seq }
                  }

            if(params.meta){
                  Channel.fromPath( "$params.meta", checkIfExists: true)
                  .set{ ch_metadata }
                  }

            preprocessing(ch_metadata, ch_seq, ch_variant)
            ch_voc=preprocessing.out.ch_voc
            ch_metadata=preprocessing.out.ch_metadata

            variant_calling(ch_voc, ch_metadata, ch_seq, ch_ref, ch_refgff, ch_reffai)
            ch_vcf=variant_calling.out.ch_vcf
            ch_stats=variant_calling.out.ch_stats
            ch_metadata=variant_calling.out.ch_metadata

            annotation(ch_vcf, ch_geneannot, ch_funcannot, ch_genecoord, ch_mutationsplit, ch_variant, ch_stats)
            ch_gvf_surveillance=annotation.out.ch_gvf_surv
            ch_variant=annotation.out.ch_variant
            ch_stats=annotation.out.ch_stats

            surveillance(ch_gvf_surveillance, ch_variant, ch_stats, ch_surveillanceIndicators, ch_metadata )

      }
      

}
