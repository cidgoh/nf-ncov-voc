#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// import modules
include { IVAR                     } from '../../modules/local/ivar'
//include { tsvTovcf             } from '../modules/local/custom'
include { BWA_INDEX                } from '../../modules/nf-core/bwa/index/main'
include { BWA_MEM                  } from '../../modules/nf-core/bwa/mem/main'
include { MINIMAP2_ALIGN           } from '../../modules/nf-core/minimap2/align/main'
include { GUNZIP                   } from '../../modules/nf-core/gunzip/main'
include { TABIX_BGZIPTABIX              } from '../../modules/nf-core/tabix/bgziptabix/main'
include { BCFTOOLS_NORM                 } from '../../modules/nf-core/bcftools/norm/main'

// import sub-workflows
include { BAM_VARIANT_CALLING_SORT_FREEBAYES_BCFTOOLS       } from '../nf-core/bam_variant_calling_sort_freebayes_bcftools/main'
include { BAM_SORT_STATS_SAMTOOLS                           } from '../nf-core/bam_sort_stats_samtools/main'


workflow VARIANT_CALLING {
    take:
        sequences_grouped
        viral_genome
        viral_genome_fai

    main:
        
        fasta = Channel.value([[params.virus_accession_id], file(viral_genome, checkIfExists: true)])  
        
        if (params.aligner == "bwa"){
            BWA_INDEX([[id: params.virus_accession_id], viral_genome])
            BWA_MEM(sequences_grouped, BWA_INDEX.out.index, false )
            bam=BWA_MEM.out.bam
        }
        else{
            MINIMAP2_ALIGN(sequences_grouped, viral_genome, true, false, false)
            bam=MINIMAP2_ALIGN.out.bam
        }
        
        BAM_SORT_STATS_SAMTOOLS(bam,fasta)
        sorted_bam=BAM_SORT_STATS_SAMTOOLS.out.bam
        bam_index=BAM_SORT_STATS_SAMTOOLS.out.bai
        bam_bai = bam_index.join(sorted_bam).map{meta, bai, bam -> [meta, bam, bai, [], [], []]}     

        if(params.ivar){
            IVAR(ch_bam.combine(ch_ref).combine(ch_refgff))
            //tsvTovcf(IVAR.out.variants)
            //ch_vcf=tsvTovcf.out.vcf

        }
        else{
            
            genome =Channel.value([[params.virus_accession_id], file(viral_genome, checkIfExists: true), file(viral_genome_fai, checkIfExists: true)])
            BAM_VARIANT_CALLING_SORT_FREEBAYES_BCFTOOLS(bam_bai , genome , [], [], [] )
            GUNZIP(BAM_VARIANT_CALLING_SORT_FREEBAYES_BCFTOOLS.out.vcf)
        }

        TABIX_BGZIPTABIX(
            GUNZIP.out.gunzip
        )

        BCFTOOLS_NORM(
            TABIX_BGZIPTABIX.out.gz_tbi, 
            viral_genome
        )
        vcf=BCFTOOLS_NORM.out.vcf
        

    emit:
        vcf
        
        

}

    