#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

script_files = "$baseDir/bin"

// include workflows
include { COVIDMVP } from './workflows/covidmvp'
include { POXMVP } from './workflows/poxmvp'
include { WASTEWATER } from './workflows/virusmvp_wastewater'

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

//if ( params.mode == 'user' && ! params.userfile ) {
//    println("When --mode user, userfile (.vcf or .fasta or .tsv) should e provided with --userfile")
//    System.exit(1)
//}


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
            if (params.infection="covid"){
                  println("Executing COVID-MVP")
                  COVIDMVP()                    
            }
            
            else (params.infection="mpox"){
                  println("Executing POX-MVP")
                  POXMVP()                      
            }
      }
      

}
