
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


process processGVCF{
  publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.vcf", mode: 'copy'
  //publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.fasta", mode: 'copy'

  tag { "ProcessGVCF" }

  input:
      path(gvcf)

  output:
      path("*.variants.vcf"), emit: vcf
      path("*.consensus.vcf")
      path("*.txt")

  script:
  """
  process_gvcf.py -d ${params.var_MinDepth} \
  -l ${params.var_MinFreqThreshold} \
  -u ${params.var_FreqThreshold} \
  -m ${gvcf.baseName}.mask.txt \
  -v ${gvcf.baseName}.variants.vcf \
  -c ${gvcf.baseName}.consensus.vcf ${gvcf}
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
        path("*.vcf"), emit: annotated_vcf

    script:
      """
      mature_peptide_annotation.py \
      --vcf_file ${peptide_vcf}\
      --annotation_file ${genome_gff}\
      --output_vcf ${peptide_vcf.baseName}.filtered.annotated.vcf
      """
}


process vcfTogvf{
  publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.gvf", mode: 'copy'
  //publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.fasta", mode: 'copy'
  tag { "VCF2GVF" }

  input:
      tuple(path(annotated_vcf), path(func_annot), path(clade_def), path(gene_coord))
      each x

  output:
      path("*gvf")

  script:
  """
    vcf2gvf.py --vcffile ${annotated_vcf}\
    --pokay ${func_annot}\
    --clades ${clade_def}\
    --gene_positions ${gene_coord}\
    --strain ${x}\
    --outvcf ${annotated_vcf.baseName}.gvf\

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
