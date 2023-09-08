process DOWNLOADPANGOALIAS{

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "scheme", mode: "copy"

    output:
        path("*.json"), emit: json
    
    script:
    """
    curl -L -o alias_key.json https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json
    """
}