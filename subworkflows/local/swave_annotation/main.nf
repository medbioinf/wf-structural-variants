include { ANNOVAR } from '../../../modules/local/annovar/main'
include { HTSLIB_BGZIPTABIX } from '../../../modules/nf-core/htslib/bgziptabix/main'
include { BCFTOOLS_FILTER as BCFTOOLS_FILTER_RARE } from '../../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_FILTER as BCFTOOLS_FILTER_INFREQUENT } from '../../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_FILTER as BCFTOOLS_FILTER_FREQUENT } from '../../../modules/nf-core/bcftools/filter/main'

workflow SWAVE_ANNOTATION {

    take:
    ch_vcf  // channel: [ meta, vcf ]

    main:
    ch_versions = channel.empty()

    ch_vcf_for_filtering = ch_vcf

    if (params.annovar_dir) {
        if (!params.annovar_db) {
            error "annovar_db must be set when annovar_dir is provided (e.g. --annovar_db hg38 or --annovar_db hs1)"
        }
        
        ch_annovar_input = ch_vcf.map { meta, vcf -> [ meta, vcf ] }

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