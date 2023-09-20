process POSTPROCESSING {

  tag {"linking output files to rsync"}

  input:
    val surveillance

  script:
      """
      rm -rf /scratch/mzanwar/COVID-MVP/nf-ncov-voc/latest_gvf/*
      ln -s ${params.outdir}/${params.prefix}/annotation_vcfTogvf/* /scratch/mzanwar/COVID-MVP/nf-ncov-voc/latest_gvf/
      
      rm -rf /scratch/mzanwar/COVID-MVP/nf-ncov-voc/latest_tsv/*
      ln -s ${params.outdir}/${params.prefix}/surveillance_surveillanceRawTsv/* /scratch/mzanwar/COVID-MVP/nf-ncov-voc/latest_tsv/
      
      rm -rf /scratch/mzanwar/COVID-MVP/nf-ncov-voc/latest_pdf/*
      ln -s ${params.outdir}/${params.prefix}/surveillance_surveillancePDF/* /scratch/mzanwar/COVID-MVP/nf-ncov-voc/latest_pdf/
      """
}