#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// import modules


include { SNPEFF_ANN                    } from '../../modules/local/snpeff_ann'
include { VCFTOGVF                      } from '../../modules/local/vcftogvf'
include { TAGPROBLEMATICSITES_NCOV      } from '../../modules/local/tagproblematicsites_ncov'
include { ANNOTATEMATPEPTIDES_NCOV      } from '../../modules/local/annotatematpeptides_ncov'

workflow ANNOTATION {
    take:
        annotation_vcf
        ch_snpeff_db
        ch_snpeff_config
        viral_genome
        ch_stats

    main:
        
        if (!params.skip_tag_problematics_sites && !params.mpox){
            problematics_sites = file(params.probvcf, checkIfExists: true)
            vcf = [ [ id:params.viral_genome_id ], [ problematics_sites ] ]
            TAGPROBLEMATICSITES_NCOV(annotation_vcf, vcf)
            annotation_vcf=TAGPROBLEMATICSITES_NCOV.out.vcf
        }        
        
        if (!params.skip_SNPEFF){
            SNPEFF_ANN (
                annotation_vcf,
                ch_snpeff_db,
                ch_snpeff_config,
                viral_genome
            )
            annotation_vcf = SNPEFF_ANN.out.vcf
        }
        

        if (!params.skip_mat_peptide_annottaion && !params.mpox){
          
            ANNOTATEMATPEPTIDES_NCOV(
              annotation_vcf,
              params.viral_gff
            )
            annotation_vcf=ANNOTATEMATPEPTIDES_NCOV.out.vcf
        }
        
        //VCF to GVF transformation
        
        json_file = file(params.genecoord, checkIfExists: true)
        json = [ [ id:params.viral_genome_id ], [ json_file ] ]
        threshold=0.75
        
        VCFTOGVF(
            annotation_vcf,
            ch_stats.map{it[1]},
            threshold,
            json, 
            true
            )
        gvf = VCFTOGVF.out.gvf

        
    emit:
        gvf
}
