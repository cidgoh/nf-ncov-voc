process articDownloadScheme{
    tag params.schemeRepoURL

    label 'internet'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "scheme", mode: "copy"

    output:
    path "${params.schemeDir}/${params.scheme}/${params.schemeVersion}/*.reference.fasta" , emit: reffasta
    path "${params.schemeDir}/${params.scheme}/${params.schemeVersion}/*.primer.bed" , emit: bed
    path "${params.schemeDir}" , emit: scheme

    script:
    """
    git clone ${params.schemeRepoURL} ${params.schemeDir}
    """
}