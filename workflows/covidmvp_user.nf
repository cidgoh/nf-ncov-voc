#!/usr/bin/env nextflow

// enable dsl2
nextflow.preview.dsl = 2
//params.ref = ".github/data/referencedb"

// import modules


include { tsvTovcf             } from '../modules/custom.nf'
include { SNPEFF               } from '../modules/snpeff.nf'
include { tagProblematicSites  } from '../modules/custom.nf'
include { annotate_mat_peptide } from '../modules/custom.nf'
include { vcfTogvf             } from '../modules/custom.nf'


workflow ncov_voc_user {
    take:
      ch_user_vcf
      ch_probvcf
      ch_geneannot
      ch_funcannot
      ch_cladedef
      ch_genecoord

    main:
      //extractMetadata(ch_metadata, ch_voc)
      //SEQKIT(extractMetadata.out.ids.combine(ch_seq))
      //BBMAP(SEQKIT.out.fasta)



      if(params.tsv){

        tsvTovcf(ch_user_vcf)
        tagProblematicSites(tsvTovcf.out.vcf.combine(ch_probvcf))
      }
      else{

        tagProblematicSites(ch_user_vcf.combine(ch_probvcf))
      }
      SNPEFF(tagProblematicSites.out.filtered_vcf)
      annotate_mat_peptide(SNPEFF.out.peptide_vcf.combine(ch_geneannot))
      vcfTogvf(annotate_mat_peptide.out.annotated_vcf.combine(ch_funcannot).combine(ch_cladedef).combine(ch_genecoord))

}
