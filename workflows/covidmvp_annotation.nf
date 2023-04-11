#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// import modules

//include { SNPEFF_BUILD                  } from '../modules/local/snpeff_build'
include { SNPEFF_ANN                    } from '../modules/local/snpeff_ann'
include { TABIX_BGZIPTABIX              } from '../modules/nf-core/tabix/bgziptabix/main'
include { BCFTOOLS_NORM                 } from '../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_INDEX                } from '../modules/nf-core/bcftools/index/main'
include { tagProblematicSites           } from '../modules/local/custom'
include { annotate_mat_peptide          } from '../modules/local/custom'
include { vcfTogvf                      } from '../modules/local/custom'


workflow annotation {
    take:
      ch_vcf
      ch_stats
      ch_snpeff_db
      ch_snpeff_config     

    main:

      tagProblematicSites(ch_vcf, "$params.prob_sites/problematic_sites_sarsCov2.vcf")
      
      TABIX_BGZIPTABIX(
        tagProblematicSites.out.filtered_vcf
      )

      
      BCFTOOLS_NORM(
        TABIX_BGZIPTABIX.out.gz_tbi, 
        params.viral_genome
      )
      snpeff_vcf=BCFTOOLS_NORM.out.vcf
      
      
      SNPEFF_ANN (
        snpeff_vcf,
        ch_snpeff_db,
        ch_snpeff_config,
        params.viral_genome
      )
      
      annotation_vcf = SNPEFF_ANN.out.vcf
      
      annotate_mat_peptide(
        annotation_vcf,
        params.viral_gff
      )

      annotated_vcf=annotate_mat_peptide.out.vcf

      
      vcfTogvf(
        annotated_vcf.combine(ch_stats, by: 0),
        params.funcannot,
        params.genecoord,
        params.mutationsplit,
        params.variant
        )

      if(params.mode == 'reference'){
        vcfTogvf.out.gvf
              .collect()
              .set{ ch_gvf_surv }
      }
      else{
        ch_gvf_surv=vcfTogvf.out.gvf
      }

    emit:
      gvf = vcfTogvf.out.gvf
      ch_stats
}
