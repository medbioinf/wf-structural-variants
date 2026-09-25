include { SWAVE_EXTRACT_VARIANT_SEQUENCES } from '../../../modules/local/swave/extract_variant_sequences/main'
include { SWAVE_SPLIT_ALLELES as SWAVE_SPLIT_VARIANT_SEQUENCES } from '../../../modules/local/swave/split_alleles/main'
include { TRF } from '../../../modules/local/trf/main'
include { SWAVE_ANNOTATE_VCF_WITH_TRF } from '../../../modules/local/swave/annotate_vcf_with_trf/main'
include { ANNOVAR } from '../../../modules/local/annovar/main'

include { HTSLIB_BGZIPTABIX } from '../../../modules/nf-core/htslib/bgziptabix/main'
include { BCFTOOLS_FILTER as BCFTOOLS_FILTER_RARE } from '../../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_FILTER as BCFTOOLS_FILTER_INFREQUENT } from '../../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_FILTER as BCFTOOLS_FILTER_FREQUENT } from '../../../modules/nf-core/bcftools/filter/main'

workflow SWAVE_ANNOTATION {

    take:
    ch_vcf          // channel: [ meta, vcf ]
    ch_gfa_fasta    // channel: [ fasta ]

    main:
    ch_versions = channel.empty()

    SWAVE_EXTRACT_VARIANT_SEQUENCES(ch_vcf, ch_gfa_fasta)
    ch_versions = ch_versions.mix(SWAVE_EXTRACT_VARIANT_SEQUENCES.out.versions_swave)

    SWAVE_SPLIT_VARIANT_SEQUENCES(SWAVE_EXTRACT_VARIANT_SEQUENCES.out.fa)
    ch_versions = ch_versions.mix(SWAVE_SPLIT_VARIANT_SEQUENCES.out.versions_swave)

    ch_split_fastas = SWAVE_SPLIT_VARIANT_SEQUENCES.out.splits.transpose()

    TRF(ch_split_fastas)
    ch_versions = ch_versions.mix(TRF.out.versions_trf)

    ch_all_dat_files = TRF.out.dat
        .map { _meta, dat -> dat }
        .collect()

    ch_all_split_fastas = ch_split_fastas
        .map { _meta, fasta -> fasta }
        .collect()

    SWAVE_ANNOTATE_VCF_WITH_TRF(ch_vcf, ch_all_dat_files, ch_all_split_fastas)
    ch_versions = ch_versions.mix(SWAVE_ANNOTATE_VCF_WITH_TRF.out.versions_swave)

    ch_vcf_for_filtering = SWAVE_ANNOTATE_VCF_WITH_TRF.out.vcf

    if (params.annovar_dir) {
        if (!params.annovar_db) {
            error "annovar_db must be set when annovar_dir is provided (e.g. --annovar_db hg38 or --annovar_db hs1)"
        }
        ch_annovar_input = ch_vcf_for_filtering.map { meta, vcf -> [ meta, vcf ] }

        ANNOVAR(
            ch_annovar_input,
            file(params.annovar_dir),
            params.annovar_db
        )
        ch_versions = ch_versions.mix(ANNOVAR.out.versions_annovar)
        ch_vcf_for_filtering = ANNOVAR.out.vcf
    }

    ch_vcf_for_filtering
        .map { meta, vcf -> [ meta, vcf, [], [] ] }
        .set { ch_htslib_input }

    HTSLIB_BGZIPTABIX(ch_htslib_input, 'compress', 'tbi', 'vcf.gz')
    ch_versions = ch_versions.mix(HTSLIB_BGZIPTABIX.out.versions_htslib)

    HTSLIB_BGZIPTABIX.out.output
        .join(HTSLIB_BGZIPTABIX.out.index)
        .set { ch_bcftools_input }

    BCFTOOLS_FILTER_RARE(ch_bcftools_input)
    ch_versions = ch_versions.mix(BCFTOOLS_FILTER_RARE.out.versions_bcftools)

    BCFTOOLS_FILTER_INFREQUENT(ch_bcftools_input)
    ch_versions = ch_versions.mix(BCFTOOLS_FILTER_INFREQUENT.out.versions_bcftools)

    BCFTOOLS_FILTER_FREQUENT(ch_bcftools_input)
    ch_versions = ch_versions.mix(BCFTOOLS_FILTER_FREQUENT.out.versions_bcftools)

    emit:
    versions = ch_versions
}