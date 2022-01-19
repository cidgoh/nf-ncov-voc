#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// import modules
include { extractVariants      } from '../modules/local/custom'
include { grabIndex            } from '../modules/local/custom'
include { BWAINDEX             } from '../modules/local/bwaindex'
include { BBMAP                } from  '../modules/local/bbmap_reformat'
include { extractMetadata      } from '../modules/local/custom'
include { SEQKIT               } from '../modules/local/seqkit'
include { SEQKITSTATS          } from '../modules/local/seqkitstats'
include { BWAMEM               } from '../modules/local/bwamem'
include { MINIMAP2             } from '../modules/local/minimap2'
include { FREEBAYES            } from '../modules/local/freebayes'
include { processGVCF          } from '../modules/local/custom'
include { BCFTOOLS             } from '../modules/local/bcftools'
include { IVAR                 } from '../modules/local/ivar'
include { tsvTovcf             } from '../modules/local/custom'




workflow variant_calling {
    take:
      ch_voc
      ch_metadata
      ch_seq
      ch_ref
      ch_refgff
      ch_reffai



    main:

      extractMetadata(ch_metadata, ch_voc)
      SEQKIT(extractMetadata.out.ids.combine(ch_seq))
      BBMAP(SEQKIT.out.fasta)
      ch_BBMap=BBMAP.out.qcfasta.unique().collect()
      SEQKITSTATS(ch_BBMap)
      ch_stats=SEQKITSTATS.out.stats

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
        ch_vcf=tsvTovcf.out.vcf

      }
      else{
        FREEBAYES(MINIMAP2.out.bam.combine(ch_ref).combine(ch_reffai),MINIMAP2.out.index)
        processGVCF(FREEBAYES.out.gvcf)
        BCFTOOLS(processGVCF.out.vcf.combine(ch_ref))
        ch_vcf=BCFTOOLS.out.normalized_vcf
      }

    emit:
      ch_vcf
      ch_stats
      ch_metadata

}
