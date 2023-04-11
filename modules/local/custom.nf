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


process IVAR_VARIANTS_TO_VCF {
    tag "$meta.id"

    conda "conda-forge::python=3.9.5 conda-forge::matplotlib=3.5.1 conda-forge::pandas=1.3.5 conda-forge::r-sys=3.4 conda-forge::regex=2021.11.10 conda-forge::scipy=1.7.3 conda-forge::biopython=1.79"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ff46c3f421ca930fcc54e67ab61c8e1bcbddfe22:1ad3da14f705eb0cdff6b5a44fea4909307524b4-0' :
        'quay.io/biocontainers/mulled-v2-ff46c3f421ca930fcc54e67ab61c8e1bcbddfe22:1ad3da14f705eb0cdff6b5a44fea4909307524b4-0' }"

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:  // This script is bundled with the pipeline, in nf-core/viralrecon/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ivar_variants_to_vcf.py \\
        $tsv \\
        ${prefix}.vcf \\
        $args \\
        > ${prefix}.variant_counts.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

process processGVCF {

  tag {"${gvcf.baseName}"}

  publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*variants.vcf", mode: 'copy'

  input:
      path(gvcf)

  
  output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml"           , emit: versions
    
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

    tag "$meta.id"

    conda "bioconda::cyvcf=0.8.0 bioconda::gffutils=0.10.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cyvcf2:0.8.0--py36_0' :
        'quay.io/biocontainers/cyvcf2:0.8.0--py36_0' }"
    
    input:
        tuple val(meta), path(vcf)
        path(prob_vcf)

    output:
        tuple val(meta), path("*.filtered.vcf"), emit: filtered_vcf

    when:
      vcf.size() > 0

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

      """
      problematic_sites_tag.py \\
        --vcffile $vcf \\
        --filter_vcf $prob_vcf \\
        --output_vcf ${prefix}.filtered.vcf \\
        $args
      """
}


process annotate_mat_peptide {

  tag "$meta.id"

  conda "bioconda::cyvcf=0.8.0 bioconda::gffutils=0.10.1"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'docker://cidgoh/nf-ncov-voc-extra:0.2' : ''}"

  input:
      tuple val(meta), path(vcf)
      path  gff

  output:
      tuple val(meta), path("*annotated.vcf"), emit: vcf
  
  when:
    vcf.size() > 0

  script:

  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mature_peptide_annotation.py \\
    --vcf_file $vcf \\
    --annotation_file $gff \\
    --output_vcf ${prefix}.annotated.vcf \\
    $args
    """
}


process vcfTogvf {

  tag "$meta.id"

  conda "bioconda::pandas=1.4.3"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.4.3' : '' }"
  
  input:
      tuple val(meta), path(vcf), path(ch_stats)
      path func_tsv
      path json
      path mutation_tsv
      path variants_tsv

  output:
      tuple val(meta), path("*.gvf"), emit: gvf

  script:

  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"

  if (params.userfile){
    input_file = file(params.userfile)
  }

  if( params.mode == 'reference'){
    """
      vcf2gvf.py --vcffile $vcf \\
      --functional_annotations $func_tsv \\ 
      --clades $variants_tsv \\
      --gene_positions $json \\
      --names_to_split $mutation_tsv \\
      --size_stats $ch_stats \\
      $args \\
      --strain ${ch_annotated_vcf.baseName.replaceAll(".qc.sorted.variants.filtered.SNPEFF.annotated", "")} --outgvf ${ch_annotated_vcf.baseName.replaceAll(".qc","_qc")}.gvf
    """
  }
  else if( params.mode == 'user' && input_file.getExtension() == "vcf"){
    """
      vcf2gvf.py --vcffile $vcf \\
      --functional_annotations $func_tsv \\
      --gene_positions $json \\
      --names_to_split $mutation_tsv \\
      $args \\
      --outgvf ${prefix}.gvf
    """
  }
  else if( params.mode == 'wastewater'){
    """
      vcf2gvf.py --vcffile $vcf \\
      --functional_annotations $func_tsv \\
      --gene_positions $json \\
      --names_to_split $mutation_tsv \\
      --size_stats $ch_stats \\
      $args \\
      --outgvf ${prefix}.gvf
    """
  }

  else{
    """
      vcf2gvf.py --vcffile $vcf \\
      --functional_annotations $func_tsv \\
      --gene_positions $json \\
      --names_to_split $mutation_tsv \\
      --size_stats $ch_stats \\
      $args \\
      --outgvf ${prefix}.gvf
    """
  }

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
