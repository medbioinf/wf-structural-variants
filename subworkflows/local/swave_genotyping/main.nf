include { SWAVE_PREDICT } from '../../../modules/local/swave/predict/main'
include { SWAVE_CALL_VARIANTS } from '../../../modules/local/swave/call_variants/main'
include { SWAVE_WRITE_VCF } from '../../../modules/local/swave/write_vcf/main'

workflow SWAVE_GENOTYPING {

    take:
    ch_dotplots     // channel: [ meta, [dotplot_files] ]
    ch_projections  // channel: [ meta, [projection_files] ]
    ch_ref_fasta    // channel: [ fasta ]
    ch_equal_paths  // channel: [ meta, txt ]

    main:
    ch_versions = channel.empty()

    SWAVE_PREDICT(ch_projections)
    ch_versions = ch_versions.mix(SWAVE_PREDICT.out.versions_swave)

    ch_calling_inputs = SWAVE_PREDICT.out.predictions
        .join(ch_projections)
        .join(ch_dotplots)
    
    SWAVE_CALL_VARIANTS(ch_calling_inputs)
    ch_versions = ch_versions.mix(SWAVE_CALL_VARIANTS.out.versions_swave)

    SWAVE_CALL_VARIANTS.out.tsv
        .map { _meta, tsv -> tsv }
        .collect()
        .map { tsv_list -> 
            def meta = [ id: 'pangenomesv' ]
            return [ meta, tsv_list ]
        }
        .set { ch_tsv_collection }
    
    SWAVE_WRITE_VCF(ch_tsv_collection, ch_ref_fasta, ch_equal_paths)
    ch_versions = ch_versions.mix(SWAVE_WRITE_VCF.out.versions_swave)

    emit:
    vcf_hap_level = SWAVE_WRITE_VCF.out.vcf.map { meta, files ->
        [ meta, files.find { f -> f.name.endsWith('hap_level.vcf') } ]
    }
    vcf_hap_level_split = SWAVE_WRITE_VCF.out.vcf.map { meta, files ->
        [ meta, files.find { f -> f.name.endsWith('hap_level.split.vcf') } ]
    }
    vcf_merged = SWAVE_WRITE_VCF.out.vcf.map { meta, files ->
        [ meta, files.find { f -> f.name.endsWith('sample_level.vcf') } ]
    }
    vcf_split = SWAVE_WRITE_VCF.out.vcf.map { meta, files ->
        [ meta, files.find { f -> f.name.endsWith('sample_level.split.vcf') } ]
    }
    versions = ch_versions
}