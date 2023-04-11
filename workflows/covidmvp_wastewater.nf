//
// Importing the modules required for the sub-workflow
//
// include { KRAKEN2_KRAKEN2               } from '../../modules/nf-core/kraken2/kraken2/main'
// include { KRAKENTOOLS_KREPORT2KRONA     } from '../../modules/nf-core/krakentools/kreport2krona/main'
include { FASTP                                         } from '../modules/nf-core/fastp/main'
include { SEQKIT_STATS                                  } from '../modules/nf-core/seqkit/stats/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_HOST         } from '../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_VIRAL        } from '../modules/nf-core/minimap2/align/main'
include { BWA_INDEX      as BWA_INDEX_HOST              } from '../modules/nf-core/bwa/index/main'
include { BWA_MEM        as BWA_MEM_HOST                } from '../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_HOST         } from '../modules/nf-core/samtools/index/main'
include { BWA_INDEX      as BWA_INDEX_VIRAL             } from '../modules/nf-core/bwa/index/main'
include { BWA_MEM        as BWA_MEM_VIRAL               } from '../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_VIRAL        } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_HOST   } from '../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_VIRAL  } from '../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_VIEW                                 } from '../modules/nf-core/samtools/view/main'
include { SAMTOOLS_FASTQ                                } from '../modules/nf-core/samtools/fastq/main'
include { SAMTOOLS_SORT                                 } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_IVAR         } from '../modules/nf-core/samtools/index/main'
include { articDownloadScheme                           } from '../modules/local/arcticdownloadscheme'
include { IVAR_TRIM                                     } from '../modules/nf-core/ivar/trim/main'
include { FREEBAYES                                     } from '../modules/nf-core/freebayes/main'
include { SAMTOOLS_MPILEUP                              } from '../modules/nf-core/samtools/mpileup/main'
include { FREYJA_UPDATE                                 } from '../modules/nf-core/freyja/update/main'
include { FREYJA_DEMIX                                  } from '../modules/nf-core/freyja/demix/main'
include { FREYJA_VARIANTS                               } from '../modules/nf-core/freyja/variants/main'
include { IVAR_VARIANTS_TO_VCF                          } from '../modules/local/custom'


workflow wastewater {
    take:
    ch_raw_reads 
    host_genome
    viral_genome


    main:
    ch_versions = Channel.empty()
    articDownloadScheme()
    if (!params.skip_fastp) {
        if (params.adapter_fasta) { adapter = file(params.adapter_fasta) 
            FASTP (
                ch_raw_reads, adapter, params.save_trimmed_fail, params.save_merged
                )
            }
        else{
            FASTP (
                ch_raw_reads, [], params.save_trimmed_fail, params.save_merged
                )  
            }
        ch_short_reads = FASTP.out.reads
        ch_versions = ch_versions.mix(FASTP.out.versions.first())
    }
    
    SEQKIT_STATS(ch_short_reads)    

    if (!params.skip_dehosting){
        if (params.aligner=='bwa') {

            BWA_INDEX_HOST(
                [[id: params.host_genome_id], host_genome]
            )
            ref_genome = BWA_INDEX_HOST.out.index
            reads = ch_short_reads
            sorted=params.bwa_sort_bam
            BWA_MEM_HOST(
                reads,
                ref_genome,
                sorted
            )
            ch_mapped_bam=BWA_MEM_HOST.out.bam
        }
        else {
            ref_genome = host_genome
            reads = ch_short_reads
            bam_format  = params.bam_format //true
            cigar_paf_format = params.cigar_paf_format //false
            cigar_bam = params.cigar_bam //false
            MINIMAP2_ALIGN_HOST ( 
                reads, 
                ref_genome, 
                bam_format, 
                cigar_paf_format, 
                cigar_bam 
            )
            ch_mapped_bam=MINIMAP2_ALIGN_HOST.out.bam
        }

        SAMTOOLS_INDEX_HOST (
            ch_mapped_bam
            )
        bam_index = SAMTOOLS_INDEX_HOST.out.bai

        SAMTOOLS_FLAGSTAT_HOST (
            ch_mapped_bam.combine (
                bam_index, by: 0
            )
            )

        SAMTOOLS_VIEW(
            ch_mapped_bam.combine(bam_index, by: 0),
            [],
            [])
        interleaved=params.interleaved
        SAMTOOLS_FASTQ(
             SAMTOOLS_VIEW.out.bam, interleaved)
        
        ch_short_reads = SAMTOOLS_FASTQ.out.fastq
    }
    
    if (params.viral_aligner=='bwa') {
        BWA_INDEX_VIRAL(
                [[id: params.viral_genome_id], viral_genome]
        )
        ref_genome = BWA_INDEX_VIRAL.out.index
        reads = ch_short_reads
        sorted=params.bwa_sort_bam
        BWA_MEM_VIRAL(
            reads,
            ref_genome,
            sorted
        )
        ch_viral_bam=BWA_MEM_VIRAL.out.bam
    }
    else {
        ref_genome = viral_genome
        reads = ch_short_reads
        bam_format  = params.bam_format //true
        cigar_paf_format = params.cigar_paf_format //false
        cigar_bam = params.cigar_bam //false
        MINIMAP2_ALIGN_VIRAL ( 
            reads, 
            ref_genome, 
            bam_format, 
            cigar_paf_format, 
            cigar_bam 
        )
        ch_viral_bam=MINIMAP2_ALIGN_VIRAL.out.bam   
    }

    SAMTOOLS_INDEX_VIRAL(
        ch_viral_bam
        )
    viral_bam_index = SAMTOOLS_INDEX_VIRAL.out.bai

    SAMTOOLS_FLAGSTAT_VIRAL(
        ch_viral_bam.combine(viral_bam_index, by: 0)
        )
    ch_ivar=ch_viral_bam.combine(viral_bam_index, by: 0)
    IVAR_TRIM(ch_ivar, articDownloadScheme.out.bed)
    SAMTOOLS_SORT(IVAR_TRIM.out.bam)
    SAMTOOLS_INDEX_IVAR(SAMTOOLS_SORT.out.bam)
    
    db_name= "freyja_db"
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

           
    FREYJA_VARIANTS(SAMTOOLS_SORT.out.bam, [[id: params.viral_genome_id], viral_genome])
    ch_freyja_variants = FREYJA_VARIANTS.out.variants
    ch_freyja_depths   = FREYJA_VARIANTS.out.depths
    FREYJA_DEMIX (
        ch_freyja_variants,
        ch_freyja_depths,
        ch_barcodes,
        ch_lineages_meta
    )
    //ch_freyja_demix = FREYJA_DEMIX.out.demix
    //ch_versions = ch_versions.mix(FREYJA_DEMIX.out.versions.first())

    IVAR_VARIANTS_TO_VCF(FREYJA_VARIANTS.out.variants)
    
    


    emit:
    stats          = SEQKIT_STATS.out.stats
    vcf            = IVAR_VARIANTS_TO_VCF.out.vcf
    variants       = FREYJA_VARIANTS.out.variants  // channel: [ val(meta), path(variants_tsv) ]
    depths         = FREYJA_VARIANTS.out.depths    // channel: [ val(meta), path(depths_tsv) ]
    demix          = FREYJA_DEMIX.out.demix        // channel: [ val(meta), path(demix_tsv) ]
    barcodes       = ch_barcodes                   // channel: [ val(meta), path(barcodes) ]
    lineages_meta  = ch_lineages_meta              // channel: [ val(meta), path(lineages_meta) ]
    //versions       = ch_versions                   // channel: [ path(versions.yml) ]
}

