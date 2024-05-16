process MERGE_LOGFILES {

    tag "$meta.id"

    conda "conda-forge::dask=2023.10.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
                'biocontainers/dask:2023.10.1-py11-ol9_cv1' : 
                'quay.io/biocontainers/dask:2023.10.1-py11-ol9_cv1' }"
    

    input:
            tuple val(meta), path(log_header)
            tuple val(meta2), path(logs)
            tuple val(meta3), path(index)
            
    output:
            tuple val(meta), path("*.tsv"), emit: merged_logfile

    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def end_date = "${params.end_date}"
    def last_index = index ? "--original_index ${index}" : ''

    """
        merge_logfiles.py \\
                --log_header ${log_header} \\
                --gvf_log ${logs} \\
                 $last_index \\
                --log_savefile ${end_date}.log.tsv
    """

}