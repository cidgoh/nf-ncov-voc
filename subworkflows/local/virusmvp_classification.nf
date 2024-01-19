nextflow.enable.dsl = 2

// Load the required modules
//
// Check input samplesheet and get read channels
//

include { PANGOLIN } from '../../modules/nf-core/pangolin'
include { MERGE_CLASSIFFICATION_METADATA } from '../../modules/local/merge_classification_report'
include { NEXTCLADE_RUN } from '../../modules/nf-core/nextclade/run'
include { NEXTCLADE_DATASETGET } from '../../modules/nf-core/nextclade/datasetget'



// Define the workflow
workflow CLASSIFICATION{
    
    take:
        metadata
        sequences

    main:
        
        // Run pangolin
        if (!params.skip_pangolin){
                PANGOLIN( sequences )
                pangolin_report = PANGOLIN.out.report
                if (metadata != 'NO_FILE'){
                    MERGE_PANGOLIN_METADATA(metadata, pangolin_report)
                    metadata = MERGE_PANGOLIN_METADATA.out.tsv
                    merged_metadata = metadata
                }
        }
        
        // RUN NEXTCLADE
        if (!params.skip_nextclade){
            if (params.virus_accession_id == "NC_045512.2"){
                dataset = 'sars-cov-2'
                reference = 'MN908947'
                tag = '2021-06-01T12:00:00Z'
            }
            else if(params.virus_accession_id == "NC_063383.1"){
                dataset = "hMPXV"
                reference = "NC_063383.1"
                tag = "2023-08-01T12:00:00Z"
            }
                NEXTCLADE_DATASETGET ( dataset, reference, tag ) 
                NEXTCLADE_RUN ( sequences, NEXTCLADE_DATASETGET.out.dataset )
                //metadata = MERGE_PANGOLIN_METADATA(metadata, PANGOLIN.out.report)
                //merged_metadata = metadata
            }
        merged_metadata = metadata

    emit:
        merged_metadata
        sequences
      
}
