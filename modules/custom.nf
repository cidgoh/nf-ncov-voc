
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

  tag { "${x}" }

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



process tsvTovcf {

    tag {"${variants_tsv.baseName.replace(".variants", "")}"}
    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.vcf", mode: 'copy'

    input:
        path(variants_tsv)

    output:
        path("*.vcf"), emit: vcf

    script:
      """
      ivar_variants_to_vcf.py ${variants_tsv} ${variants_tsv.baseName}.vcf
      """
}

process tagProblematicSites {

    tag {"${vcf.baseName.replace(".variants", "")}"}
    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.vcf", mode: 'copy'

    input:
        tuple(path(vcf), path(prob_vcf))

    output:
        path("*.vcf"), emit: filtered_vcf

    script:
      """
      problematic_sites_tag.py \
      --vcffile ${vcf} \
      --filter_vcf ${prob_vcf} \
      --output_vcf ${vcf.baseName}.filtered.vcf
      """
}

process annotate_mat_peptide {

    tag {"${peptide_vcf.baseName.replace(".variants.filtered.annotated", "")}"}
    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.vcf", mode: 'copy'

    input:
        tuple(path(peptide_vcf), path(genome_gff))

    output:
        path("*.vcf")

    script:
      """
      mature_peptide_annotation.py \
      --vcf_file ${peptide_vcf}\
      --annotation_file ${genome_gff}\
      --output_vcf ${peptide_vcf.baseName}.filtered.vcf
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



process vcfTotsv {

    tag {"vcfTotsv${annotated_vcf.baseName}"}
    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.tsv", mode: 'copy'

    input:
        path(annotated_vcf)

    output:
        path("*.tsv")

    script:
      """
      vcf2tsv.py ${annotated_vcf} ${annotated_vcf.baseName}.tsv
      """
}

process merge_vcfs {
  publishDir path: "${params.outdir}/freebayes"
  //cpus params.cpus.toInteger()

  input:
  file(vcf) from vcf_ch.collect()
  file(vcfidx) from vcfidx_ch.collect()
  file(fasta)
  file(faidx)
  file(vcftxt) from vcftxt_ch

  output:
  file("${params.project}.vcf.gz")
  file("${params.project}.vcf.gz.tbi")

  script:
  """
  gunzip -cd \$(cat $vcftxt) | vcffirstheader | bgzip -c > ${params.project}_dirty.vcf.gz
  gsort ${params.project}_dirty.vcf.gz $faidx | vcfuniq \
      | bgzip -c > ${params.project}_dirty_sorted.vcf.gz
  bcftools norm -c all -f $fasta --multiallelics - --threads ${task.cpus} \
      --output ${params.project}.vcf.gz --output-type z \
      ${params.project}_dirty_sorted.vcf.gz
  tabix -p vcf ${params.project}.vcf.gz
  """
}
