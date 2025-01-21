#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// import modules


include { SNPEFF_ANN               } from '../../modules/local/snpeff_ann'
include { VCFTOGVF                 } from '../../modules/local/vcftogvf'
include { TAGPROBLEMATICSITES_NCOV } from '../../modules/local/tagproblematicsites_ncov'
include { ANNOTATEMATPEPTIDES_NCOV } from '../../modules/local/annotatematpeptides_ncov'

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
    wastewater_data = []

    ch_viral_genome = Channel
        .fromPath(viral_genome)
        .ifEmpty { error("Cannot find viral genome file: ${viral_genome}") }
        .map { file -> [[id: params.virus_accession_id], file] }
        .collect()

    if (!params.skip_problematics_sites && params.virus_accession_id == "NC_045512.2") {
        problematics_sites = file(params.probvcf, checkIfExists: true)
        vcf = [[id: params.virus_accession_id], [problematics_sites]]
        TAGPROBLEMATICSITES_NCOV(annotation_vcf, vcf)
        annotation_vcf = TAGPROBLEMATICSITES_NCOV.out.vcf
    }
    ch_viral_genome_value = ch_viral_genome.first()
    if (!params.skip_SNPEFF) {
        SNPEFF_ANN(
            annotation_vcf,
            ch_snpeff_db.collect(),
            ch_snpeff_config.collect(),
            ch_viral_genome_value,
        )
        annotation_vcf = SNPEFF_ANN.out.vcf
    }

    ch_json_value = ch_json.first()
    if (!params.skip_peptide_annottaion && params.virus_accession_id == "NC_045512.2") {

        ANNOTATEMATPEPTIDES_NCOV(
            annotation_vcf,
            ch_json_value,
        )
        annotation_vcf = ANNOTATEMATPEPTIDES_NCOV.out.vcf
    }

    //VCF to GVF transformation

    /*if (ch_json == []){
            json_file = file(params.genecoord, checkIfExists: true)
            json = [ [ id:params.virus_accession_id ], [ json_file ] ]
        }
        else{
            json = ch_json
        }*/
    threshold = 0.75
    if (params.wastewater) {
        lineage = []
        wastewater_data = true
        data_description = "Wastewater"

        // For wastewater, combine annotation_vcf with its corresponding stats file
        combined_input = annotation_vcf.join(ch_stats, by: 0)
    }
    else {
        data_description = "Clinical"
        wastewater_data = false
        lineage = true

        if (ch_stats) {
            ch_stats_value = ch_stats.collect().map { it -> it[1] }.first()
            combined_input = annotation_vcf.combine(ch_stats_value) 
        }
        else {
            combined_input = annotation_vcf.map { meta, vcf ->
                [meta, vcf, []]
            }
        }
    }

    ch_json_value = ch_json.first()

    VCFTOGVF(
        combined_input,
        threshold,
        ch_json_value,
        lineage,
        wastewater_data,
        data_description,
    )

    gvf = VCFTOGVF.out.gvf

    emit:
    gvf
}
