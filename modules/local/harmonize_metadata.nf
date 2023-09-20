process METADATA_HARMONIZER {
    tag "$meta.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://pandas/pandas:pip-all': '' }"

    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.gz", mode: 'copy'

    input:
        tuple val(meta), path(metadata)
        path(yaml)
        val(dataset)
        
    output:
        tuple val(meta), path("*.csv.gz"), emit: gz

    script:

        """
        harmonize_metadata.py \\
        --metadata $metadata \\
        --config $yaml \\
        --data $dataset \\
        --outfile harmonized_metadata.csv.gz

        """
}    
