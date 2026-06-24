include { SWAVE_PREDICT } from '../../../modules/local/swave_predict/main'

workflow SWAVE_GENOTYPING {

    take:
    ch_projections

    main:
    ch_versions = channel.empty()

    SWAVE_PREDICT(ch_projections)
    ch_versions = ch_versions.mix(SWAVE_PREDICT.out.versions)

    emit:
    predictions = SWAVE_PREDICT.out.predictions
    versions = ch_versions
}