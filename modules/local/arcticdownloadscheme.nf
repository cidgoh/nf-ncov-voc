process articDownloadScheme {
    tag params.schemeRepoURL

    label 'internet'

    output:
    path "${params.schemeDir}/${params.scheme}/${params.schemeVersion}/*.reference.fasta", emit: reffasta
    path "${params.schemeDir}/${params.scheme}/${params.schemeVersion}/*.primer.bed", emit: bed
    path "${params.schemeDir}", emit: scheme

    script:
    """
    git clone ${params.schemeRepoURL} ${params.schemeDir}
    """
}
