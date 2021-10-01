#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// import modules

include { grabIndex            } from '../modules/custom.nf'
include { BWAINDEX             } from '../modules/bwaindex.nf'
include { BBMAP                } from  '../modules/bbmap_reformat.nf'
include { extractMetadata      } from '../modules/custom.nf'
include { SEQKIT               } from '../modules/seqkit.nf'
include { BWAMEM               } from '../modules/bwamem.nf'
include { MINIMAP2             } from '../modules/minimap2.nf'
include { FREEBAYES            } from '../modules/freebayes.nf'
include { processGVCF          } from '../modules/custom.nf'
include { BCFTOOLS             } from '../modules/bcftools.nf'
include { IVAR                 } from '../modules/ivar.nf'
include { tsvTovcf             } from '../modules/custom.nf'
include { SNPEFF               } from '../modules/snpeff.nf'
include { tagProblematicSites  } from '../modules/custom.nf'
include { annotate_mat_peptide } from '../modules/custom.nf'
include { vcfTogvf             } from '../modules/custom.nf'


workflow ncov_voc {
    take:
      ch_seq
      ch_metadata
      ch_voc
      ch_ref
      ch_refgff
      ch_reffai
      ch_probvcf
      ch_geneannot
      ch_funcannot
      ch_cladedef
      ch_genecoord

    main:
      if(!params.single_genome){
        extractMetadata(ch_metadata, ch_voc)
        SEQKIT(extractMetadata.out.ids.combine(ch_seq))

        BBMAP(SEQKIT.out.fasta)
      }
      else{
        BBMAP(ch_seq)
      }

      if (params.bwa){
        if ( params.bwa_index ){
        grabIndex("${params.bwa_index}")
        grabIndex.out
                .set{ ch_index }
      } else {
        BWAINDEX(ch_ref)
        BWAINDEX.out
                .set{ ch_index }
              }
        BWAMEM(BBMAP.out.qcfasta.combine(ch_ref), ch_index)
      }
      else{
        MINIMAP2(BBMAP.out.qcfasta.combine(ch_ref))
      }
      if(params.ivar){
        IVAR(MINIMAP2.out.bam.combine(ch_ref).combine(ch_refgff))
        tsvTovcf(IVAR.out.variants)
        tagProblematicSites(tsvTovcf.out.vcf.combine(ch_probvcf))
      }
      else{
        FREEBAYES(MINIMAP2.out.bam.combine(ch_ref).combine(ch_reffai),MINIMAP2.out.index)
        processGVCF(FREEBAYES.out.gvcf)
        BCFTOOLS(processGVCF.out.vcf.combine(ch_ref))
        tagProblematicSites(BCFTOOLS.out.normalized_vcf.combine(ch_probvcf))
      }
      SNPEFF(tagProblematicSites.out.filtered_vcf)
      annotate_mat_peptide(SNPEFF.out.peptide_vcf.combine(ch_geneannot))
      vcfTogvf(annotate_mat_peptide.out.annotated_vcf.combine(ch_funcannot).combine(ch_cladedef).combine(ch_genecoord))

}
