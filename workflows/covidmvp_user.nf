#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// import modules


include { tsvTovcf             } from '../modules/custom.nf'
include { MINIMAP2             } from '../modules/minimap2.nf'
include { FREEBAYES            } from '../modules/freebayes.nf'
include { BCFTOOLS             } from '../modules/bcftools.nf'
include { processGVCF          } from '../modules/custom.nf'
include { SNPEFF               } from '../modules/snpeff.nf'
include { tagProblematicSites  } from '../modules/custom.nf'
include { annotate_mat_peptide } from '../modules/custom.nf'
include { vcfTogvf             } from '../modules/custom.nf'


workflow ncov_voc_user {
    take:
      ch_probvcf
      ch_geneannot
      ch_funcannot
      ch_cladedef
      ch_genecoord

    main:
      //extractMetadata(ch_metadata, ch_voc)
      //SEQKIT(extractMetadata.out.ids.combine(ch_seq))
      //BBMAP(SEQKIT.out.fasta)


      if(params.userfile){
        Channel.fromPath( "$params.userfile", checkIfExists: true)
             .set{ ch_user_file }
      }

      if(params.input_type == 'fasta'){

        Channel.fromPath( "$params.refdb/MN908947.3.fasta", checkIfExists: true)
              .set{ ch_ref }

        Channel.fromPath( "$params.refdb/*.fai", checkIfExists: true)
              .set{ ch_reffai }

        MINIMAP2(ch_user_file.combine(ch_ref))
        FREEBAYES(MINIMAP2.out.bam.combine(ch_ref).combine(ch_reffai),MINIMAP2.out.index)
        processGVCF(FREEBAYES.out.gvcf)
        BCFTOOLS(processGVCF.out.vcf.combine(ch_ref))
        tagProblematicSites(BCFTOOLS.out.normalized_vcf.combine(ch_probvcf))
      }

      else if(params.input_type == 'tsv'){

        tsvTovcf(ch_user_file)
        tagProblematicSites(tsvTovcf.out.vcf.combine(ch_probvcf))
      }

      else{
        Channel.fromPath( "$params.userfile", checkIfExists: true)
              .set{ ch_user_file }

        tagProblematicSites(ch_user_file.combine(ch_probvcf))
      }

      SNPEFF(tagProblematicSites.out.filtered_vcf)
      annotate_mat_peptide(SNPEFF.out.peptide_vcf.combine(ch_geneannot))
      vcfTogvf(annotate_mat_peptide.out.annotated_vcf.combine(ch_funcannot).combine(ch_cladedef).combine(ch_genecoord))

}
