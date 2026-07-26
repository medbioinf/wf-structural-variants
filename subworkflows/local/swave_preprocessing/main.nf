include { SWAVE_EXTRACT_ALLELES } from '../../../modules/local/swave/extract_alleles/main'
include { SWAVE_SPLIT_ALLELES } from '../../../modules/local/swave/split_alleles/main'
include { SWAVE_GENERATE_DOTPLOTS } from '../../../modules/local/swave/generate_dotplots/main'
include { SWAVE_GENERATE_PROJECTIONS } from '../../../modules/local/swave/generate_projections/main'

workflow SWAVE_PREPROCESSING {

    take:
    ch_bed
    ch_gfa_fasta
    ch_ref_fasta

    main:
    ch_versions = channel.empty()

    SWAVE_EXTRACT_ALLELES(ch_bed, ch_gfa_fasta.toList())
    ch_versions = ch_versions.mix(SWAVE_EXTRACT_ALLELES.out.versions)

    ch_sample_alleles = SWAVE_EXTRACT_ALLELES.out.fa
        .filter { meta, _fa -> meta.is_ref != true }
    
    SWAVE_SPLIT_ALLELES(ch_sample_alleles)
    ch_versions = ch_versions.mix(SWAVE_SPLIT_ALLELES.out.versions)

    ch_dotplot_inputs = SWAVE_SPLIT_ALLELES.out.splits
        .transpose()
        .map { meta, fa ->
            def new_meta = meta.clone()
            new_meta.sample = meta.id
            new_meta.id = fa.baseName
            return [ new_meta, fa ]
        }

    SWAVE_GENERATE_DOTPLOTS(ch_dotplot_inputs, ch_ref_fasta.toList(), ch_gfa_fasta.toList())
    ch_versions = ch_versions.mix(SWAVE_GENERATE_DOTPLOTS.out.versions)

    SWAVE_GENERATE_PROJECTIONS(SWAVE_GENERATE_DOTPLOTS.out.dotplots)
    ch_versions = ch_versions.mix(SWAVE_GENERATE_PROJECTIONS.out.versions)

    emit:
    alleles_fasta = SWAVE_SPLIT_ALLELES.out.splits
    dotplots = SWAVE_GENERATE_DOTPLOTS.out.dotplots
    dotplot_pngs = SWAVE_GENERATE_DOTPLOTS.out.pngs
    projections = SWAVE_GENERATE_PROJECTIONS.out.projections
    projections_pngs = SWAVE_GENERATE_PROJECTIONS.out.pngs
    versions = ch_versions
}