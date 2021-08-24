#!/usr/bin/env nextflow

// enable dsl2
nextflow.preview.dsl = 2

params.seq = "$baseDir/.github/data/sequence"
params.meta = "$baseDir/.github/data/metadata"
params.voc = "$baseDir/.github/data/VOC"
params.refdb = "$baseDir/.github/data/refdb"
params.ref_gff = "$baseDir/.github/data/features"
params.prob_sites = "$baseDir/.github/data/problematic_sites"
params.genome_annotation = "$baseDir/.github/data/genome_annotation"
params.functional_annotation = "$baseDir/.github/data/functional_annotation"
params.clade_defining_mutations = "$baseDir/.github/data/clade_defining"
params.gene_coordinates = "$baseDir/.github/data/gene_coordinates"



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

      //Channel.from('B.1.1.7', 'B.1.351', 'B.1.351.2', 'B.1.351.3', 'P.1', 'P.1.1', 'P.1.2', 'P.1.4', 'P.1.6', 'P.1.7', 'B.1.617.2', 'AY.1', 'AY.2', 'AY.3', 'AY.3.1',  'B.1.525', 'B.1.526', 'B.1.617.1', 'C.37')
      //      .set{ch_voc}
      //Channel.from('B.1.351.3', 'B.1.351.2', 'B.1.351', 'AY.3.1', 'AY.3', 'AY.2', 'AY.1', 'B.1.617.2')
      //      .set{ch_voc}
      Channel.from('B.1.351.3')
            .set{ ch_voc }

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

      Channel.fromPath( "$params.prob_sites/*.vcf", checkIfExists: true)
            .set{ ch_probvcf }

      Channel.fromPath( "$params.genome_annotation/*.gff", checkIfExists: true)
            .set{ ch_geneannot }

      Channel.fromPath( "$params.functional_annotation/*.tsv", checkIfExists: true)
            .set{ ch_funcannot }

      Channel.fromPath( "$params.clade_defining_mutations/*.tsv", checkIfExists: true)
            .set{ ch_cladedef }

      Channel.fromPath( "$params.gene_coordinates/*.json", checkIfExists: true)
            .set{ ch_genecoord }


   main:

       //println("This will call Illumina workflow")
       ncov_voc(ch_seq, ch_metadata, ch_voc, ch_ref, ch_refgff, ch_reffai, ch_probvcf, ch_geneannot, ch_funcannot, ch_cladedef, ch_genecoord)
}
