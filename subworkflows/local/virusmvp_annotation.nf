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
        ch_json

    main:
        lineage = true
        wastewater = []
        
        if (!params.skip_tag_problematics_sites && params.virus_accession_id == "NC_045512.2"){
            problematics_sites = file(params.probvcf, checkIfExists: true)
            vcf = [ [ id:params.virus_accession_id ], [ problematics_sites ] ]
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
        

        if (!params.skip_mat_peptide_annottaion && params.virus_accession_id == "NC_045512.2"){
          
            ANNOTATEMATPEPTIDES_NCOV(
              annotation_vcf,
              ch_json
            )
            annotation_vcf=ANNOTATEMATPEPTIDES_NCOV.out.vcf
        }
        
        //VCF to GVF transformation
        
        if (ch_json == []){
            json_file = file(params.genecoord, checkIfExists: true)
            json = [ [ id:params.virus_accession_id ], [ json_file ] ]
        }
        else{
            json = ch_json
        }
        threshold=0.75
        if (params.wastewater){
            lineage = []
            wastewater = true
            data_description = "Wastewater"
        }
        else{
            data_description = "Clinical"
        }

        VCFTOGVF(
            annotation_vcf,
            ch_stats.map{it[1]},
            threshold,
            json, 
            lineage, 
            wastewater,
            data_description,
            )
        gvf = VCFTOGVF.out.gvf

        
    emit:
        gvf
}
