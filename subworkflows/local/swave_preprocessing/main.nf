include { SWAVE_EXTRACT_ALLELES } from '../../../modules/local/swave_extract_alleles/main'
include { SEQKIT_SORT } from '../../../modules/nf-core/seqkit/sort/main'
include { SEQKIT_SPLIT2 as SEQKIT_SPLIT_BY_SIZE } from '../../../modules/nf-core/seqkit/split2/main'
include { SEQKIT_SPLIT2 as SEQKIT_SPLIT_BY_LENGTH } from '../../../modules/nf-core/seqkit/split2/main'
include { SWAVE_GENERATE_DOTPLOTS } from '../../../modules/local/swave_generate_dotplots/main'
include { SWAVE_GENERATE_PROJECTIONS } from '../../../modules/local/swave_generate_projections/main'

workflow SWAVE_PREPROCESSING {

    take:
    ch_bed
    ch_gfa_fasta
    ch_ref_fasta

    main:
    ch_versions = channel.empty()
    ch_fasta_for_split = channel.empty()
    ch_dotplot_inputs = channel.empty()

    SWAVE_EXTRACT_ALLELES(ch_bed, ch_gfa_fasta.toList())
    ch_versions = ch_versions.mix(SWAVE_EXTRACT_ALLELES.out.versions)

    SEQKIT_SORT(SWAVE_EXTRACT_ALLELES.out.fa)
    ch_versions = ch_versions.mix(SEQKIT_SORT.out.versions_seqkit)

    ch_fasta_for_split = SEQKIT_SORT.out.fastx
    .map { meta, fa ->
        def new_meta = meta.clone()
        new_meta.single_end = true  // single-end for compatibility with seqkit split
        return [ new_meta, fa ]
    }

    SEQKIT_SPLIT_BY_SIZE(ch_fasta_for_split)
    ch_versions = ch_versions.mix(SEQKIT_SPLIT_BY_SIZE.out.versions_seqkit)

    ch_intermediate_splits = SEQKIT_SPLIT_BY_SIZE.out.reads
        .transpose()
        .map { meta, fa ->
            def new_meta = meta.clone()
            new_meta.sample = meta.id
            new_meta.id = fa.baseName
            return [ new_meta, fa ]
        }
    
    SEQKIT_SPLIT_BY_LENGTH(ch_intermediate_splits)
    ch_versions = ch_versions.mix(SEQKIT_SPLIT_BY_LENGTH.out.versions_seqkit)

    ch_dotplot_inputs = SEQKIT_SPLIT_BY_LENGTH.out.reads
        .map { meta, files -> [ meta, files ] }
        .transpose()
        .map { meta, fa ->
            def new_meta = meta.clone()
            new_meta.sample = meta.sample ?: meta.id
            new_meta.id = fa.baseName
            return [ new_meta, fa ]
        }

    SWAVE_GENERATE_DOTPLOTS(ch_dotplot_inputs, ch_ref_fasta.toList(), ch_gfa_fasta.toList())
    ch_versions = ch_versions.mix(SWAVE_GENERATE_DOTPLOTS.out.versions)

    SWAVE_GENERATE_PROJECTIONS(SWAVE_GENERATE_DOTPLOTS.out.dotplots)
    ch_versions = ch_versions.mix(SWAVE_GENERATE_PROJECTIONS.out.versions)

    emit:
    alleles_fasta = SWAVE_EXTRACT_ALLELES.out.fa
    dotplots = SWAVE_GENERATE_DOTPLOTS.out.dotplots
    dotplot_pngs = SWAVE_GENERATE_DOTPLOTS.out.pngs
    projections = SWAVE_GENERATE_PROJECTIONS.out.projections
    projections_pngs = SWAVE_GENERATE_PROJECTIONS.out.pngs
    versions = ch_versions
}