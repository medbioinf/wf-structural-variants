include { SWAVE_PREDICT } from '../../../modules/local/swave_predict/main'
include { SWAVE_CALL_VARIANTS } from '../../../modules/local/swave_call_variants/main'
include { SWAVE_WRITE_VCF } from '../../../modules/local/swave_write_vcf/main'

workflow SWAVE_GENOTYPING {

    take:
    ch_dotplots
    ch_projections
    ch_reference_fasta

    main:
    ch_versions = channel.empty()

    SWAVE_PREDICT(ch_projections)
    ch_versions = ch_versions.mix(SWAVE_PREDICT.out.versions)

    ch_calling_inputs = SWAVE_PREDICT.out.predictions
        .join(ch_projections)
        .join(ch_dotplots)
        
    
    SWAVE_CALL_VARIANTS(ch_calling_inputs)
    ch_versions = ch_versions.mix(SWAVE_CALL_VARIANTS.out.versions)

    SWAVE_CALL_VARIANTS.out.tsv
        .map { _meta, tsv -> tsv }
        .collect()
        .map { tsv_list -> 
            def meta = [ id: 'swave_population' ]
            return [ meta, tsv_list ]
        }
        .set { ch_tsv_collection }
    
    SWAVE_WRITE_VCF(ch_tsv_collection, ch_reference_fasta)
    ch_versions = ch_versions.mix(SWAVE_WRITE_VCF.out.versions)

    emit:
    predictions = SWAVE_PREDICT.out.predictions
    variants_tsv = SWAVE_CALL_VARIANTS.out.tsv
    vcf = SWAVE_WRITE_VCF.out.vcf
    versions = ch_versions
}