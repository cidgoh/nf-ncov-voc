#!/usr/bin/env nextflow

// enable dsl2
nextflow.preview.dsl = 2

params.seq = ".github/data/sequence"
params.meta = ".github/data/metadata"
params.voc = ".github/data/VOC"
params.refdb = ".github/data/refdb"
params.ref_gff = ".github/data/features"

// include modules
include {printHelp} from './modules/help.nf'
//include {makeFastqSearchPath} from './modules/util.nf'

// import subworkflows
include {ncov_voc} from './workflows/covidmvp.nf'


if (params.help){
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

      Channel.from('B.1.1.7','B.1.351', 'B.1.351.2', 'B.1.351.3', 'P.1', 'P.1.1', 'P.1.2', 'B.1.617.2', 'AY.1', 'AY.2', 'B.1.525', 'B.1.526', 'B.1.617.1', 'C.37')
            .set{ch_voc}
      //Channel.from('B.1.351', 'B.1.525')
      //      .set{ch_voc}

      Channel.fromPath( "$params.seq/*.fasta", checkIfExists: true)
	         .set{ ch_seq }

      Channel.fromPath( "$params.meta/metadata*.tsv", checkIfExists: true)
            .set{ ch_metadata }

      Channel.fromPath( "$params.ref_gff/*.gff3", checkIfExists: true)
            .set{ ch_refgff }

      Channel.fromPath( "$params.refdb/*.fai", checkIfExists: true)
            .set{ ch_reffai }

      Channel.fromPath( "$params.refdb/MN908947.3.fasta", checkIfExists: true)
            .set{ ch_ref }

   main:

       //println("This will call Illumina workflow")
       ncov_voc(ch_seq, ch_metadata, ch_voc, ch_ref, ch_refgff, ch_reffai)
}
