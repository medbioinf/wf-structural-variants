include { SWAVE_PREDICT } from '../../../modules/local/swave_predict/main'
include { SWAVE_CALL_VARIANTS } from '../../../modules/local/swave_call_variants/main'

workflow SWAVE_GENOTYPING {

    take:
    ch_dotplots
    ch_projections

    main:
    ch_versions = channel.empty()

    SWAVE_PREDICT(ch_projections)
    ch_versions = ch_versions.mix(SWAVE_PREDICT.out.versions)

    ch_calling_inputs = SWAVE_PREDICT.out.predictions
        .join(ch_projections)
        .join(ch_dotplots)
        
    
    SWAVE_CALL_VARIANTS(ch_calling_inputs)
    ch_versions = ch_versions.mix(SWAVE_CALL_VARIANTS.out.versions)

    emit:
    predictions = SWAVE_PREDICT.out.predictions
    variants_tsv = SWAVE_CALL_VARIANTS.out.tsv
    versions = ch_versions
}