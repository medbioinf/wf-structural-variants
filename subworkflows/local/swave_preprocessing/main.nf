include { SWAVE_EXTRACT_ALLELES } from '../../../modules/local/swave_extract_alleles/main'

workflow SWAVE_PREPROCESSING {

    take:
    ch_bed
    ch_gfa_fasta

    main:
    ch_versions = channel.empty()

    SWAVE_EXTRACT_ALLELES(ch_bed, ch_gfa_fasta)
    ch_versions = ch_versions.mix(SWAVE_EXTRACT_ALLELES.out.versions)

    emit:
    alleles_fa = SWAVE_EXTRACT_ALLELES.out.fa
    versions = ch_versions

}