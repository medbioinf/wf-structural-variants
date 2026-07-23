include { HTSLIB_BGZIPTABIX } from '../../../modules/nf-core/htslib/bgziptabix/main'                                                                                                        
include { BCFTOOLS_FILTER as BCFTOOLS_FILTER_RARE } from '../../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_FILTER as BCFTOOLS_FILTER_INFREQUENT } from '../../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_FILTER as BCFTOOLS_FILTER_FREQUENT } from '../../../modules/nf-core/bcftools/filter/main'

workflow SWAVE_ANNOTATION {

    take:
    ch_vcf_split

    main:
    ch_versions = channel.empty()

    ch_vcf_split
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
    vcf_rare = BCFTOOLS_FILTER_RARE.out.vcf
    vcf_infrequent = BCFTOOLS_FILTER_INFREQUENT.out.vcf
    vcf_frequent = BCFTOOLS_FILTER_FREQUENT.out.vcf
    versions = ch_versions
}
