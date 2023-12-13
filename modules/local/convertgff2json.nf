/*
    This module takes a GFF file as input and converts it to JSON format using the gff2json script.
*/

// Define the process to run
process CONVERTGFFTOJSON {
     tag "$meta.id"

    conda "bioconda::pandas=1.4.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/pandas:1.4.3' : ''}"
    
    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.json", mode: 'copy'

    input:
        tuple val(meta), path(gff)
        path(color)
        path(alias)
    
    output:
        tuple val(meta), path("*.json"), emit: json

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        def color = color ? "--color_file ${color}" : ''
        def alias = alias ? "--alias_file ${alias}" : ''
        
        """
        gff2json.py \\
        --gff_file ${gff} \\
        $color \\
        $alias \\
        --json_file ${prefix}.json
        
        """
}
