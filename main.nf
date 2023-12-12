#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

script_files = "$baseDir/bin"

// include workflows
include { COVIDMVP } from './workflows/covidmvp'
//include { POXMVP } from './workflows/poxmvp'
include { WASTEWATER } from './workflows/virusmvp_wastewater'
//include { FLUMVP                    } from './workflows/flumvp'

// include subworkflows
include { CONFIGURE_VIRUSMVP } from './subworkflows/local/virusmvp_config'

// include modules
include {printHelp              } from './modules/local/help'
include {cidgohHeader           } from './modules/local/header'
include {workflowHeader         } from './modules/local/wf_header'



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
      if (params.wastewater){
            WASTEWATER()
      }
      else{
             // Building configuration file for COVID-MVP
            if(!params.skip_configuration){
                  gff_file = file(params.viral_gff, checkIfExists: true)
                  gff = [ [ id:params.virus_accession_id ], gff_file ]
                  CONFIGURE_VIRUSMVP(gff)
                  ch_json=CONFIGURE_VIRUSMVP.out.ch_json
                  ch_snpeff_db=CONFIGURE_VIRUSMVP.out.ch_snpeff_db
                  ch_snpeff_config=CONFIGURE_VIRUSMVP.out.ch_snpeff_config
                  ch_viral_fai=CONFIGURE_VIRUSMVP.out.ch_viral_fai
            }
            if (params.virus_accession_id = "NC_045512.2"){
                  println("Executing COVID-MVP")
                  COVIDMVP(ch_json, ch_snpeff_db, ch_snpeff_config, ch_viral_fai) 
                  
            }

            else if (params.virus_accession_id = "NC_063383.1"){
                  println("Executing POX-MVP")
                  //POXMVP(ch_json, ch_snpeff_db, ch_snpeff_config, ch_viral_fai)                    
            }
            else {
                  println("Executing FLU-MVP")
                  //FLUMVP()                                              
            }
      }
      

}
