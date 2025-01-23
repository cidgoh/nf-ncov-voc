#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// import modules
include { IVAR                                        } from '../../modules/local/ivar'
include { IVAR_VARIANTS_TO_VCF                        } from '../../modules/local/archive/custom'
include { BWA_INDEX                                   } from '../../modules/nf-core/bwa/index/main'
include { BWA_MEM                                     } from '../../modules/nf-core/bwa/mem/main'
include { MINIMAP2_ALIGN                              } from '../../modules/nf-core/minimap2/align/main'
include { GUNZIP                                      } from '../../modules/nf-core/gunzip/main'
include { TABIX_BGZIPTABIX                            } from '../../modules/nf-core/tabix/bgziptabix/main'
include { BCFTOOLS_NORM                               } from '../../modules/nf-core/bcftools/norm/main'

// import sub-workflows
include { BAM_VARIANT_CALLING_SORT_FREEBAYES_BCFTOOLS } from '../nf-core/bam_variant_calling_sort_freebayes_bcftools/main'
include { BAM_SORT_STATS_SAMTOOLS                     } from '../nf-core/bam_sort_stats_samtools/main'


workflow VARIANT_CALLING {
    take:
    sequences_grouped
    viral_genome
    viral_genome_fai

    main:

    // Create a channel for the viral genome file
    ch_viral_genome = Channel
        .fromPath(viral_genome)
        .ifEmpty { error("Cannot find viral genome file: ${viral_genome}") }
        .map { file -> [[id: params.virus_accession_id], file] }
        .collect()

    ch_viral_genome_fai = Channel
        .fromPath(viral_genome_fai)
        .ifEmpty { error("Cannot find viral genome file: ${viral_genome}") }
        .collect()


    if (params.viral_aligner == "bwa") {
        BWA_INDEX(ch_viral_genome)
        index = BWA_INDEX.out.index
        BWA_MEM(sequences_grouped, index, ch_viral_genome.map { it -> it[1] }, true)
        ch_bam = BWA_MEM.out.bam
    }
    else {
        MINIMAP2_ALIGN(sequences_grouped, ch_viral_genome, true, false, false)
        ch_bam = MINIMAP2_ALIGN.out.bam
    }

    // Modify BAM_SORT_STATS_SAMTOOLS call to match expected inputs
    BAM_SORT_STATS_SAMTOOLS(ch_bam, ch_viral_genome)
    sorted_bam = BAM_SORT_STATS_SAMTOOLS.out.bam
    bam_index = BAM_SORT_STATS_SAMTOOLS.out.bai

    bam_bai = sorted_bam
        .join(bam_index)
        .map { meta, bam, bai ->
            tuple(meta, bam, bai)
        }

    if (params.ivar) {
        IVAR(ch_bam.combine(ch_viral_genome).combine(params.ch_refgff))
        IVAR_VARIANTS_TO_VCF(IVAR.out.variants)
        vcf = IVAR_VARIANTS_TO_VCF.out.vcf
    }
    else {
        ch_input = bam_bai.map { meta, bam, bai ->
            tuple(meta, bam, bai, [], [], [])
        }

        ch_viral_genome_with_index = ch_viral_genome
            .combine(ch_viral_genome_fai)
            .map { meta, fasta, fai ->
                [meta[0], fasta, fai[0]]
            }

        // Create a value channel from ch_viral_genome_with_index
        ch_viral_genome_with_index_value = ch_viral_genome_with_index.first()
        BAM_VARIANT_CALLING_SORT_FREEBAYES_BCFTOOLS(
            ch_input,
            ch_viral_genome_with_index_value,
            [[], []],
            [[], []],
            [[], []]
        )
        GUNZIP(BAM_VARIANT_CALLING_SORT_FREEBAYES_BCFTOOLS.out.vcf)
        vcf = GUNZIP.out.gunzip
    }

    TABIX_BGZIPTABIX(vcf)

    BCFTOOLS_NORM(
        TABIX_BGZIPTABIX.out.gz_tbi,
        ch_viral_genome,
    )
    vcf = BCFTOOLS_NORM.out.vcf

    emit:
    vcf
}
