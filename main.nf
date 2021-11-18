#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// path to required files
params.refdb = "$baseDir/.github/data/refdb"
params.ref_gff = "$baseDir/.github/data/features"
params.prob_sites = "$baseDir/.github/data/problematic_sites"
params.genome_annotation = "$baseDir/.github/data/genome_annotation"
params.functional_annotation = "$baseDir/.github/data/functional_annotation"
//params.clade_defining_mutations = "$baseDir/.github/data/clade_defining"
params.gene_coordinates = "$baseDir/.github/data/gene_coordinates"
params.mutation_names = "$baseDir/.github/data/multi_aa_names"


// include modules
include {printHelp              } from './modules/help.nf'
include {cidgohHeader           } from './modules/header.nf'
include {workflowHeader         } from './modules/wf_header.nf'
include {extractVariants        } from './modules/custom.nf'
include { virusseqMapLineage    } from './modules/custom.nf'

// import subworkflows
include {ncov_voc         } from './workflows/covidmvp.nf'
include {ncov_voc_user    } from './workflows/covidmvp_user.nf'

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
if ( params.data == 'virusseq' && ! params.virusseq_meta ) {
    println("VirusSeq metadata file is required when data type is virusSeq, provide a .tsv file with --virusseq")
    System.exit(1)
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
      println(params.mode+ " activated!!\n\n\n")

      if(params.variants){
        Channel.fromPath( "$params.variants", checkIfExists: true)
              .set{ ch_variants }
          }
      else{
        Channel.fromPath( "$baseDir/.github/data/variants/*.tsv", checkIfExists: true)
              .set{ ch_variant }
          }

      if(params.mode == 'user'){

        Channel.fromPath( "$params.prob_sites/*.vcf", checkIfExists: true)
              .set{ ch_probvcf }

        Channel.fromPath( "$params.genome_annotation/*.gff", checkIfExists: true)
              .set{ ch_geneannot }

        Channel.fromPath( "$params.functional_annotation/*.tsv", checkIfExists: true)
              .set{ ch_funcannot }

        Channel.fromPath( "$params.gene_coordinates/*.json", checkIfExists: true)
              .set{ ch_genecoord }

        Channel.fromPath( "$params.mutation_names/*.tsv", checkIfExists: true)
              .set{ ch_mutationsplit }


        ncov_voc_user( ch_probvcf, ch_geneannot, ch_funcannot, ch_variant, ch_genecoord, ch_mutationsplit)
      }


      else if(params.mode == 'reference'){

        if(params.seq){
          Channel.fromPath( "$params.seq", checkIfExists: true)
    	         .set{ ch_seq }
        }
        else{
          Channel.fromPath( "$baseDir/.github/data/sequence/*.fasta", checkIfExists: true)
    	         .set{ ch_seq }
        }

        if(params.meta){
          Channel.fromPath( "$params.meta", checkIfExists: true)
               .set{ ch_metadata }
        }

        else{
          Channel.fromPath( "$baseDir/.github/data/metadata/*.tsv", checkIfExists: true)
                .set{ ch_metadata }
        }

        if(params.data == "virusseq"){
          if(params.virusseq_meta){

          Channel.fromPath( "$params.virusseq_meta", checkIfExists: true)
               .set{ ch_virusseq_metadata }
             }

          virusseqMapLineage(ch_virusseq_metadata.combine(ch_metadata))
          virusseqMapLineage.out.mapped
              .set{ch_metadata}
          //extractVariants(ch_variant.combine(virusseqMapLineage.out.mapped))
          //extractVariants.out.lineages
          //    .splitText()
          //    .set{ ch_voc }

        }

        extractVariants(ch_variant.combine(ch_metadata))
        extractVariants.out.lineages
                .splitText()
                .set{ ch_voc }


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

        Channel.fromPath( "$params.gene_coordinates/*.json", checkIfExists: true)
              .set{ ch_genecoord }

        Channel.fromPath( "$params.mutation_names/*.tsv", checkIfExists: true)
              .set{ ch_mutationsplit }


        ncov_voc(ch_voc, ch_metadata, ch_seq, ch_ref, ch_refgff, ch_reffai, ch_probvcf, ch_geneannot, ch_funcannot, ch_genecoord, ch_mutationsplit, ch_variant )
      }

}
