
process grabIndex {

    tag { "Grab_Covid-19_index" }

    input:
      path(index_folder)

    output:
      file("*.fasta.*")

    script:
      """
      ln -sf $index_folder/*.fasta.* ./
      """
}



process extractMetadata {
  publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.tsv", mode: 'copy'
  //publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.fasta", mode: 'copy'

  tag { "ExtractMetadata_${x}" }

  input:
      path(metadata)
      //path(sequence)
      each x

  output:
      path("*.tsv")
      path("*.txt"), emit: ids

  script:
    """
    extract_VOC_GISAID.py \
    --table ${metadata} \
    --voc ${x} \
    --startdate ${params.startdate} \
    --enddate ${params.enddate}
    """
}

process processGVCF{
  publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.vcf", mode: 'copy'
  //publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.fasta", mode: 'copy'

  tag { "ProcessGVCF" }

  input:
      path(variants)

  output:
      path("*.variants.vcf"), emit: vcf
      path("*.consensus.vcf")
      path("*.txt")

  script:
  """
  process_gvcf.py -d ${params.var_MinDepth} \
                        -l ${params.var_MinFreqThreshold} \
                        -u ${params.var_FreqThreshold} \
                        -m ${variants.baseName}.mask.txt \
                        -v ${variants.baseName}.variants.vcf \
                        -c ${variants.baseName}.consensus.vcf ${variants}
  """
}
