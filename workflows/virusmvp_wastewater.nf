// Importing the modules required for the sub-workflow

include { FASTP                 as WW_FASTP                         } from '../modules/nf-core/fastp/main'
include { SEQKIT_STATS          as WW_SEQKIT_STATS                  } from '../modules/nf-core/seqkit/stats/main'
include { MINIMAP2_ALIGN        as WW_MINIMAP2_ALIGN_HOST           } from '../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN        as WW_MINIMAP2_ALIGN_VIRAL          } from '../modules/nf-core/minimap2/align/main'
include { BWA_INDEX             as WW_BWA_INDEX_HOST                } from '../modules/nf-core/bwa/index/main'
include { BWA_INDEX             as WW_BWA_INDEX_VIRAL               } from '../modules/nf-core/bwa/index/main'
include { SAMTOOLS_INDEX        as WW_SAMTOOLS_INDEX_HOST           } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX        as WW_SAMTOOLS_INDEX_VIRAL          } from '../modules/nf-core/samtools/index/main'
include { BWA_MEM               as WW_BWA_MEM_HOST                  } from '../modules/nf-core/bwa/mem/main'
include { BWA_MEM               as WW_BWA_MEM_VIRAL                 } from '../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_FLAGSTAT     as WW_SAMTOOLS_FLAGSTAT_HOST        } from '../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_FLAGSTAT     as WW_SAMTOOLS_FLAGSTAT_VIRAL       } from '../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_VIEW         as WW_SAMTOOLS_VIEW                 } from '../modules/nf-core/samtools/view/main'
include { SAMTOOLS_FASTQ        as WW_SAMTOOLS_FASTQ                } from '../modules/nf-core/samtools/fastq/main'
include { SAMTOOLS_SORT         as WW_SAMTOOLS_SORT                 } from '../modules/nf-core/samtools/sort/main'
//include { SAMTOOLS_INDEX        as WW_SAMTOOLS_INDEX_IVAR              } from '../../modules/nf-core/samtools/index/main'
include { articDownloadScheme                                       } from '../modules/local/arcticdownloadscheme'
include { IVAR_TRIM             as WW_IVAR_TRIM                     } from '../modules/nf-core/ivar/trim/main'
include { FREEBAYES             as WW_FREEBAYES                     } from '../modules/nf-core/freebayes/main'
include { SNPEFF_BUILD                                              } from '../modules/local/snpeff_build'
include { SAMTOOLS_FAIDX                                            } from '../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_MPILEUP      as WW_SAMTOOLS_MPILEUP              } from '../modules/nf-core/samtools/mpileup/main'
include { FREYJA_UPDATE                                             } from '../modules/nf-core/freyja/update/main'
include { FREYJA_DEMIX                                              } from '../modules/nf-core/freyja/demix/main'
include { FREYJA_VARIANTS                                           } from '../modules/nf-core/freyja/variants/main'
include { FREYJA_BOOT                                           } from '../modules/nf-core/freyja/boot/main'
include { IVAR_VARIANTS_TO_VCF  as WW_IVAR_VARIANTS_TO_VCF          } from '../modules/local/custom'

include {INPUT_CHECK            } from '../subworkflows/local/input_check'
include { BAM_VARIANT_DEMIX_BOOT_FREYJA } from '../subworkflows/nf-core/bam_variant_demix_boot_freyja/main'


workflow WASTEWATER {
    
    main:

        if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
                ch_versions = Channel.empty()
                INPUT_CHECK (
                    ch_input
                )
        ch_raw_reads=INPUT_CHECK.out.reads
    // Building SNPEFF database for COVID
        ch_snpeff_db     = Channel.empty()
        ch_snpeff_config = Channel.empty()
        ch_voc = Channel.empty()
        
        SNPEFF_BUILD (
                params.viral_genome,
                params.viral_gbk
        )
        ch_snpeff_db     = SNPEFF_BUILD.out.db
        ch_snpeff_config = SNPEFF_BUILD.out.config

        SAMTOOLS_FAIDX([[id: params.viral_genome_id], params.viral_genome], [[],[]])
        ch_viral_fai=SAMTOOLS_FAIDX.out.fai
    
        ch_versions = Channel.empty()
        articDownloadScheme()
        if (!params.skip_fastp) {
            if (params.adapter_fasta) { adapter = file(params.adapter_fasta) 
                WW_FASTP (
                    ch_raw_reads, adapter, params.save_trimmed_fail, params.save_merged
                    )
                }
            else{
                WW_FASTP (
                    ch_raw_reads, [], params.save_trimmed_fail, params.save_merged
                    )  
                }
            ch_short_reads = WW_FASTP.out.reads
            ch_versions = ch_versions.mix(WW_FASTP.out.versions.first())
        }
        
        WW_SEQKIT_STATS(ch_short_reads)    

        if (!params.skip_dehosting){
            if (params.aligner=='bwa') {

                WW_BWA_INDEX_HOST(
                    [[id: params.host_genome_id], params.host_genome]
                )
                ref_genome = WW_BWA_INDEX_HOST.out.index
                reads = ch_short_reads
                sorted=params.bwa_sort_bam
                WW_BWA_MEM_HOST(
                    reads,
                    ref_genome,
                    sorted
                )
                ch_mapped_bam=WW_BWA_MEM_HOST.out.bam
            }
            else {
                ref_genome = host_genome
                reads = ch_short_reads
                bam_format  = params.bam_format //true
                cigar_paf_format = params.cigar_paf_format //false
                cigar_bam = params.cigar_bam //false
                WW_MINIMAP2_ALIGN_HOST ( 
                    reads, 
                    ref_genome, 
                    bam_format, 
                    cigar_paf_format, 
                    cigar_bam 
                )
                ch_mapped_bam=WW_MINIMAP2_ALIGN_HOST.out.bam
            }

            WW_SAMTOOLS_INDEX_HOST (
                ch_mapped_bam
                )
            bam_index = WW_SAMTOOLS_INDEX_HOST.out.bai

            WW_SAMTOOLS_FLAGSTAT_HOST (
                ch_mapped_bam.combine (
                    bam_index, by: 0
                )
                )

            WW_SAMTOOLS_VIEW(
                ch_mapped_bam.combine(bam_index, by: 0),
                [],
                [])
            interleaved=params.interleaved
            WW_SAMTOOLS_FASTQ(
                WW_SAMTOOLS_VIEW.out.bam, interleaved)
            
            ch_short_reads = WW_SAMTOOLS_FASTQ.out.fastq
        }
        
        if (params.viral_aligner=='bwa') {
            WW_BWA_INDEX_VIRAL(
                    [[id: params.viral_genome_id], params.viral_genome]
            )
            ref_genome = WW_BWA_INDEX_VIRAL.out.index
            reads = ch_short_reads
            sorted=params.bwa_sort_bam
            WW_BWA_MEM_VIRAL(
                reads,
                ref_genome,
                sorted
            )
            ch_viral_bam=WW_BWA_MEM_VIRAL.out.bam
        }
        else {
            ref_genome = params.viral_genome
            reads = ch_short_reads
            bam_format  = params.bam_format //true
            cigar_paf_format = params.cigar_paf_format //false
            cigar_bam = params.cigar_bam //false
            WW_MINIMAP2_ALIGN_VIRAL ( 
                reads, 
                ref_genome, 
                bam_format, 
                cigar_paf_format, 
                cigar_bam 
            )
            ch_viral_bam=WW_MINIMAP2_ALIGN_VIRAL.out.bam   
        }

        WW_SAMTOOLS_INDEX_VIRAL(
            ch_viral_bam
            )
        viral_bam_index = WW_SAMTOOLS_INDEX_VIRAL.out.bai

        WW_SAMTOOLS_FLAGSTAT_VIRAL(
            ch_viral_bam.combine(viral_bam_index, by: 0)
            )
        ch_ivar=ch_viral_bam.combine(viral_bam_index, by: 0)
        WW_IVAR_TRIM(ch_ivar, articDownloadScheme.out.bed)
        WW_SAMTOOLS_SORT(WW_IVAR_TRIM.out.bam)
        //WW_SAMTOOLS_INDEX_VIRAL(WW_SAMTOOLS_SORT.out.bam)
        
        db_name= "freyja_db"
        //BAM_VARIANT_DEMIX_BOOT_FREYJA(WW_SAMTOOLS_SORT.out.bam, [[id: params.viral_genome_id], params.viral_genome], 100, db_name)

        //ch_stats=wastewater.out.stats        
        //annotation(wastewater.out.vcf, ch_stats, ch_snpeff_db, ch_snpeff_config)
        
        FREYJA_UPDATE(db_name)
        FREYJA_UPDATE
                .out
                .barcodes
                .map { barcodes  -> [ [], barcodes ] }
                .set { ch_barcodes }

        FREYJA_UPDATE
            .out
            .lineages_meta
            .map { lineages  -> [ [], lineages ] }
            .set { ch_lineages_meta }

            
        FREYJA_VARIANTS(WW_SAMTOOLS_SORT.out.bam, [[id: params.viral_genome_id], params.viral_genome])
        ch_freyja_variants = FREYJA_VARIANTS.out.variants
        ch_freyja_depths   = FREYJA_VARIANTS.out.depths
        FREYJA_DEMIX (
            ch_freyja_variants,
            ch_freyja_depths,
            ch_barcodes,
            ch_lineages_meta
        )
        ch_freyja_demix = FREYJA_DEMIX.out.demix
        ch_versions = ch_versions.mix(FREYJA_DEMIX.out.versions.first())
        
        WW_IVAR_VARIANTS_TO_VCF(FREYJA_VARIANTS.out.variants)
    
    


    emit:
        stats          = WW_SEQKIT_STATS.out.stats
        vcf            = WW_IVAR_VARIANTS_TO_VCF.out.vcf
        variants       = FREYJA_VARIANTS.out.variants  // channel: [ val(meta), path(variants_tsv) ]
        //depths         = FREYJA_VARIANTS.out.depths    // channel: [ val(meta), path(depths_tsv) ]
        demix          = FREYJA_DEMIX.out.demix        // channel: [ val(meta), path(demix_tsv) ]
        barcodes       = ch_barcodes                   // channel: [ val(meta), path(barcodes) ]
        lineages_meta  = ch_lineages_meta              // channel: [ val(meta), path(lineages_meta) ]
        versions       = ch_versions                   // channel: [ path(versions.yml) ]
}

