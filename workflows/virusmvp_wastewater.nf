// Importing the modules required for the sub-workflow

include { FASTP as WW_FASTP                               } from '../modules/nf-core/fastp/main'
include { SEQKIT_STATS as WW_SEQKIT_STATS                 } from '../modules/nf-core/seqkit/stats/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_HOST           } from '../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_VIRAL          } from '../modules/nf-core/minimap2/align/main'
include { BWA_INDEX as WW_BWA_INDEX_HOST                  } from '../modules/nf-core/bwa/index/main'
include { BWA_INDEX as WW_BWA_INDEX_VIRAL                 } from '../modules/nf-core/bwa/index/main'
include { SAMTOOLS_INDEX as WW_SAMTOOLS_INDEX_HOST        } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as WW_SAMTOOLS_INDEX_VIRAL       } from '../modules/nf-core/samtools/index/main'
include { BWA_MEM as WW_BWA_MEM_HOST                      } from '../modules/nf-core/bwa/mem/main'
include { BWA_MEM as WW_BWA_MEM_VIRAL                     } from '../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_FLAGSTAT as WW_SAMTOOLS_FLAGSTAT_HOST  } from '../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_FLAGSTAT as WW_SAMTOOLS_FLAGSTAT_VIRAL } from '../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_VIEW as WW_SAMTOOLS_VIEW               } from '../modules/nf-core/samtools/view/main'
include { SAMTOOLS_FASTQ as WW_SAMTOOLS_FASTQ             } from '../modules/nf-core/samtools/fastq/main'
include { SAMTOOLS_SORT as WW_SAMTOOLS_SORT               } from '../modules/nf-core/samtools/sort/main'
//include { SAMTOOLS_INDEX        as WW_SAMTOOLS_INDEX_IVAR              } from '../../modules/nf-core/samtools/index/main'
include { articDownloadScheme                             } from '../modules/local/arcticdownloadscheme'
include { IVAR_TRIM as WW_IVAR_TRIM                       } from '../modules/nf-core/ivar/trim/main'
include { FREEBAYES as WW_FREEBAYES                       } from '../modules/nf-core/freebayes/main'
include { SNPEFF_BUILD                                    } from '../modules/local/snpeff_build'
include { SAMTOOLS_FAIDX                                  } from '../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_MPILEUP as WW_SAMTOOLS_MPILEUP         } from '../modules/nf-core/samtools/mpileup/main'
include { FREYJA_UPDATE                                   } from '../modules/nf-core/freyja/update/main'
include { FREYJA_DEMIX                                    } from '../modules/nf-core/freyja/demix/main'
include { FREYJA_VARIANTS                                 } from '../modules/nf-core/freyja/variants/main'
include { IVAR_VARIANTS_TO_VCF as WW_IVAR_VARIANTS_TO_VCF } from '../modules/local/archive/custom'
include { XZ_DECOMPRESS                                   } from '../modules/nf-core/xz/decompress/main'
include { SURVEILLANCE                                    } from '../subworkflows/local/virusmvp_surveillance'
include { INPUT_CHECK                                     } from '../subworkflows/local/input_check'
include { BAM_VARIANT_DEMIX_BOOT_FREYJA                   } from '../subworkflows/nf-core/bam_variant_demix_boot_freyja/main'
include { ANNOTATION                                      } from '../subworkflows/local/virusmvp_annotation'
include { GVF_PROCESSING_ANNOTATION                       } from '../subworkflows/local/virusmvp_gvf_processing_annotations'

workflow WASTEWATER {
    take:
    ch_json
    ch_snpeff_db
    ch_snpeff_config

    main:
    ch_short_reads = Channel.empty()
    metadata = Channel.empty()
    if (params.input) {
        ch_input = file(params.input)
    }
    else {
        exit(1, 'Input samplesheet not specified!')
    }
    ch_versions = Channel.empty()
    INPUT_CHECK(
        ch_input
    )
    ch_short_reads = INPUT_CHECK.out.reads


    //ch_voc = Channel.empty()
    ch_versions = Channel.empty()
    articDownloadScheme()
    if (!params.skip_fastp) {
        if (params.adapter_fasta) {
            adapter = file(params.adapter_fasta)
            WW_FASTP(
                ch_short_reads,
                adapter,
                params.save_trimmed_fail,
                params.save_merged,
            )
        }
        else {
            WW_FASTP(
                ch_short_reads,
                [],
                params.save_trimmed_fail,
                params.save_merged,
            )
        }
        ch_short_reads = WW_FASTP.out.reads
        ch_versions = ch_versions.mix(WW_FASTP.out.versions)
    }

    WW_SEQKIT_STATS(ch_short_reads)

    viral_genome = [
        [id: params.viral_genome_id, single_end: true],
        file(params.viral_genome, checkIfExists: true),
    ]
    if (!params.skip_dehosting) {
        human_genome = [
            [id: params.host_genome_id, single_end: true],
            file(params.host_genome, checkIfExists: true),
        ]
        if (params.dehost_aligner == 'bwa') {
            // Check if BWA index files already exist
            if (!bwaIndexExists(params.host_genome)) {
                WW_BWA_INDEX_HOST(human_genome)
                human_genome_index = WW_BWA_INDEX_HOST.out.index
            }
            else {
                // If index files exist, create a channel with the existing index
                human_genome_index = Channel
                    .fromPath("${params.host_genome}.*")
                    .collect()
                    .map { files ->
                        [human_genome[0], files]
                    }
            }
            reads = ch_short_reads
            sorted = params.bwa_sort_bam
            WW_BWA_MEM_HOST(
                reads,
                human_genome_index,
                human_genome,
                sorted,
            )
            ch_mapped_bam = WW_BWA_MEM_HOST.out.bam
        }
        else {

            reads = ch_short_reads
            bam_format = params.bam_format
            //true
            cigar_paf_format = params.cigar_paf_format
            //false
            cigar_bam = params.cigar_bam
            //false
            MINIMAP2_ALIGN_HOST(
                reads,
                human_genome,
                bam_format,
                cigar_paf_format,
                cigar_bam,
            )
            ch_mapped_bam = MINIMAP2_ALIGN_HOST.out.bam
        }

        WW_SAMTOOLS_INDEX_HOST(
            ch_mapped_bam
        )
        bam_index = WW_SAMTOOLS_INDEX_HOST.out.bai

        WW_SAMTOOLS_FLAGSTAT_HOST(
            ch_mapped_bam.combine(
                bam_index,
                by: 0
            )
        )

        WW_SAMTOOLS_VIEW(
            ch_mapped_bam.combine(bam_index, by: 0),
            [[], []],
            [],
        )
        interleaved = params.interleaved
        WW_SAMTOOLS_FASTQ(
            WW_SAMTOOLS_VIEW.out.bam,
            interleaved,
        )

        ch_short_reads = WW_SAMTOOLS_FASTQ.out.fastq
    }

    if (params.viral_aligner == 'bwa') {
        WW_BWA_INDEX_VIRAL(viral_genome)
        viral_genome_index = WW_BWA_INDEX_VIRAL.out.index
        reads = ch_short_reads
        sorted = params.bwa_sort_bam
        WW_BWA_MEM_VIRAL(
            reads,
            viral_genome_index,
            viral_genome,
            sorted,
        )
        ch_viral_bam = WW_BWA_MEM_VIRAL.out.bam
    }
    else {

        reads = ch_short_reads
        bam_format = params.bam_format
        //true
        cigar_paf_format = params.cigar_paf_format
        //false
        cigar_bam = params.cigar_bam
        //false
        MINIMAP2_ALIGN_VIRAL(
            reads,
            viral_genome,
            bam_format,
            cigar_paf_format,
            cigar_bam,
        )
        ch_viral_bam = MINIMAP2_ALIGN_VIRAL.out.bam
    }

    WW_SAMTOOLS_INDEX_VIRAL(
        ch_viral_bam
    )
    viral_bam_index = WW_SAMTOOLS_INDEX_VIRAL.out.bai

    WW_SAMTOOLS_FLAGSTAT_VIRAL(
        ch_viral_bam.combine(viral_bam_index, by: 0)
    )
    ch_ivar = ch_viral_bam.combine(viral_bam_index, by: 0)
    WW_IVAR_TRIM(ch_ivar, articDownloadScheme.out.bed)
    WW_SAMTOOLS_SORT(WW_IVAR_TRIM.out.bam, viral_genome)


    FREYJA_UPDATE(params.db_name)
    ch_barcodes = FREYJA_UPDATE.out.barcodes
    ch_lineages_meta = FREYJA_UPDATE.out.lineages_meta

    BAM_VARIANT_DEMIX_BOOT_FREYJA(WW_SAMTOOLS_SORT.out.bam, file(params.viral_genome), 100, params.db_name, ch_barcodes, ch_lineages_meta)
    WW_IVAR_VARIANTS_TO_VCF(BAM_VARIANT_DEMIX_BOOT_FREYJA.out.variants)
    vcf = WW_IVAR_VARIANTS_TO_VCF.out.vcf
    ch_stats = WW_SEQKIT_STATS.out.stats

    ANNOTATION(vcf, ch_snpeff_db, ch_snpeff_config, params.viral_genome, ch_stats, ch_json)
    annotatted_vcf = ANNOTATION.out.gvf

    GVF_PROCESSING_ANNOTATION(annotatted_vcf)
    annotated_gvf = GVF_PROCESSING_ANNOTATION.out.annotation_gvf
    if (!params.skip_surveillance) {
        SURVEILLANCE(annotated_gvf, metadata)
    }

    emit:
    versions = ch_versions // channel: [ path(versions.yml) ]
}

def bwaIndexExists(genome_file) {
    def index_suffixes = ['.amb', '.ann', '.bwt', '.pac', '.sa']
    return index_suffixes.every { suffix ->
        file("${genome_file}${suffix}").exists()
    }
}
