process mergePangolinMetadata{

    tag {"Merge Pangolin report into Metadata"}

    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.tsv", mode: 'copy'

    input:
      tuple (path(metadata), path(pangolin_report))

    output:
      path("Metadata_lineage.tsv"), emit: lineage_assigned

    script:
      """
      merge_pangolin_metadata.py --metadata ${metadata} --pangolin ${pangolin_report} --output Metadata_lineage.tsv

      """
}

process virusseqMapLineage{

    tag {"Mapping VirusSeq dataset to GISAID for lineages"}

    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.tsv", mode: 'copy'

    input:
      tuple (path(metadata), path(gisaid_metadata))

    output:
      path("VirusSeq_Mapped.tsv"), emit: mapped

    script:
      """
      map_virusseq_GISAID.py --virusseq ${metadata} --gisaid ${gisaid_metadata} --output VirusSeq_Mapped.tsv

      """
}

process extractVariants {
  tag {"Extracting VOCs, VOIs and VUMs"}
  input:
      tuple (path(variants), path(metadata))
  output:
      path("*.txt"), emit: lineages

  script:

    """
    parse_variants.py \
    --variants ${variants} \
    --metadata ${metadata} \
    --outfile metada_lineages.txt
    """
}

process grabIndex {

    tag { "grabing_SARS-CoV-2_index" }

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
    tag { "Extracting Metadata and IDS for VOCs, VOIs, & VUMs" }

    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.tsv.gz", mode: 'copy'

    input:
      path(metadata)
      each x

    output:
      path("*.tsv.gz")
      path("*.txt"), emit: ids

    script:
      """
      extract_metadata.py --startdate ${params.startdate} --enddate ${params.enddate} --table ${metadata} --voc ${x}
      """
}


process tsvTovcf {

    tag {"${variants_tsv.baseName}"}

    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.vcf", mode: 'copy'

    input:
        path(variants_tsv)

    output:
        path("*.vcf"), emit: vcf

    script:
      """
      ivar_variants_to_vcf.py ${variants_tsv} ${variants_tsv}.vcf
      """
}


process processGVCF {

  tag {"${gvcf.baseName}"}

  publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*variants.vcf", mode: 'copy'

  input:
      path(gvcf)

  output:
      path("*.variants.vcf"), emit: vcf
      path("*.consensus.vcf")
      path("*.txt")

  when:
      gvcf.size()>0

  script:
      """
      process_gvcf.py -d ${params.var_MinDepth} \
      -l ${params.lower_ambiguityFrequency} \
      -u ${params.upper_ambiguityFrequency} \
      -m ${gvcf.baseName}.mask.txt \
      -v ${gvcf.baseName}.variants.vcf \
      -c ${gvcf.baseName}.consensus.vcf ${gvcf}
      """
}


process tagProblematicSites {
    tag {"${vcf.baseName}"}

    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.vcf", mode: 'copy'

    input:
        tuple(path(vcf), path(prob_vcf))

    output:
        path("*.vcf"), emit: filtered_vcf

    when:
      vcf.size() > 0

    script:
      """
      problematic_sites_tag.py \
      --vcffile ${vcf} \
      --filter_vcf ${prob_vcf} \
      --output_vcf ${vcf.baseName}.filtered.vcf
      """
}


process annotate_mat_peptide {

    tag {"${peptide_vcf.baseName}"}

    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.vcf", mode: 'copy'

    input:
        tuple(path(peptide_vcf), path(genome_gff))

    output:
        path("*.vcf"), emit: annotated_vcf

    when:
      peptide_vcf.size() > 0

    script:
      """
      mature_peptide_annotation.py \
      --vcf_file ${peptide_vcf}\
      --annotation_file ${genome_gff}\
      --output_vcf ${peptide_vcf.baseName}.annotated.vcf
      """
}


process vcfTogvf {

  tag {"${ch_annotated_vcf.baseName}"}

  publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.gvf", mode: 'copy'

  input:
      tuple(path(ch_annotated_vcf), path(func_annot), path(gene_coord), path(mutation_split), path(variants), path(stats))

  output:
      path("*.gvf"), emit: gvf


  script:
  if (params.userfile){
    input_file = file(params.userfile)
  }

  if( params.mode == 'reference'){
    """
      vcf2gvf.py --vcffile ${ch_annotated_vcf} \
      --functional_annotations ${func_annot}  \
      --clades ${variants} \
      --gene_positions ${gene_coord} \
      --names_to_split ${mutation_split} \
      --size_stats ${stats} \
      --strain ${ch_annotated_vcf.baseName.replaceAll(".qc.sorted.variants.filtered.SNPEFF.annotated", "")} --outgvf ${ch_annotated_vcf.baseName.replaceAll(".qc","_qc")}.gvf
    """
  }
  else if( params.mode == 'user' && input_file.getExtension() == "vcf"){
    """
      vcf2gvf.py --vcffile ${ch_annotated_vcf} \
      --functional_annotations ${func_annot}  \
      --gene_positions ${gene_coord} \
      --names_to_split ${mutation_split} \
      --outgvf ${ch_annotated_vcf.baseName}.gvf
    """
  }

  else{
    """
      vcf2gvf.py --vcffile ${ch_annotated_vcf} \
      --functional_annotations ${func_annot}  \
      --gene_positions ${gene_coord} \
      --names_to_split ${mutation_split} \
      --size_stats ${stats} \
      --outgvf ${ch_annotated_vcf.baseName}.gvf
    """
  }



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

process surveillanceRawTsv {

  tag {"Generating Raw Surveillance Data (TSV)"}

  publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.tsv", mode: 'copy'

  input:
      path(gvf)
      tuple(path(variants), path(stats))

  output:
      path ('*.tsv'), emit: surveillancetsv

  script:
  if( params.mode == 'user' ){
    """
    gvf2tsv.py --gvf_files ${gvf} \
    --user
    """
  }
  else{
    """
    gvf2tsv.py --gvf_files ${gvf} \
    --clades ${variants} \
    --table ${stats} \
    --all_variants

    """
  }


}

process surveillancePDF {

  tag {"Generating Surveillance reports PDF"}

  publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.pdf", mode: 'copy'

  input:
      each tsv
      path(surveillanceindicators)
      path(metadata)

  output:
      path("*.pdf"), emit: surveillance_pdf

  script:
  if( params.mode == 'reference' ){

    """
    surveillance_report_pdf.py --tsv ${tsv} \
    --functions_table ${surveillanceindicators} \
    --metadata ${metadata} \
    --virusseq True > ${tsv.baseName}.tex

    tectonic -X compile ${tsv.baseName}.tex --reruns 3 --keep-intermediates --keep-logs


    """
  }
  else{
    """
    surveillance_report_pdf.py --tsv ${tsv} \
    --functions_table ${surveillanceindicators} > ${tsv.baseName}.tex 
    
    tectonic -X compile ${tsv.baseName}.tex --reruns 3 --keep-intermediates --keep-logs

    """
  }
}
