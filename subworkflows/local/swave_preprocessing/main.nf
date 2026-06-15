include { SWAVE_EXTRACT_ALLELES } from '../../../modules/local/swave_extract_alleles/main'
include { SWAVE_GENERATE_DOTPLOTS } from '../../../modules/local/swave_generate_dotplots/main'

workflow SWAVE_PREPROCESSING {

    take:
    ch_bed
    ch_gfa_fasta
    ch_ref_fasta

    main:
    ch_versions = channel.empty()

    SWAVE_EXTRACT_ALLELES(ch_bed, ch_gfa_fasta.toList())
    ch_versions = ch_versions.mix(SWAVE_EXTRACT_ALLELES.out.versions)

    SWAVE_GENERATE_DOTPLOTS(SWAVE_EXTRACT_ALLELES.out.fa, ch_ref_fasta.toList(), ch_gfa_fasta.toList())
    ch_versions = ch_versions.mix(SWAVE_GENERATE_DOTPLOTS.out.versions)

    emit:
    alleles_fasta = SWAVE_EXTRACT_ALLELES.out.fa
    dotplot_matrices = SWAVE_GENERATE_DOTPLOTS.out.dotplot_matrices
    dotplot_pngs = SWAVE_GENERATE_DOTPLOTS.out.pngs
    versions = ch_versions
}