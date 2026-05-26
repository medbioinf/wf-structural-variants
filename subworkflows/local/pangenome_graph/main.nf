include { MINIGRAPH_CONSTRUCT } from '../../../modules/local/minigraph_construct/main'
include { MINIGRAPH_CALL } from '../../../modules/local/minigraph_call/main'
include { GFATOOLS_GFA2FA } from '../../../modules/local/gfatools_gfa2fa/main'

workflow PANGENOME_GRAPH {

    take:
    ch_reference
    ch_assemblies

    main:
    ch_versions = channel.empty()
    ch_gfa_with_meta = channel.empty()
    ch_info = channel.empty()
    ch_fa = channel.empty()

    // TODO: Potentially add alternative graph construction/calling methods here (e.g. Cactus, PGGB, etc.) and allow user to select via params

    if (!params.gfa) {

        ch_assemblies
            .multiMap { meta, fasta ->
                metas: meta
                fastas: fasta
            }
            .set { ch_split_assemblies }

        ch_minigraph_assemblies = ch_split_assemblies.fastas.collect()
        MINIGRAPH_CONSTRUCT(ch_reference, ch_minigraph_assemblies)

        ch_gfa_with_meta = MINIGRAPH_CONSTRUCT.out.gfa
        ch_info = MINIGRAPH_CONSTRUCT.out.info
        ch_versions = ch_versions.mix(MINIGRAPH_CONSTRUCT.out.versions)
    } else {
        ch_gfa_with_meta = channel.fromPath(params.gfa).map { gfa -> [ [id: gfa.baseName], gfa] }
    }

    if (!params.gfa) {
        GFATOOLS_GFA2FA(ch_gfa_with_meta)
        ch_fa = GFATOOLS_GFA2FA.out.fa.map { _meta, fa -> fa }
        ch_versions = ch_versions.mix(GFATOOLS_GFA2FA.out.versions)
    } else {
        ch_fa = channel.fromPath(params.gfa2fa_fa)
    }

    ch_gfa_raw = ch_gfa_with_meta.map { _meta, gfa -> gfa }
    MINIGRAPH_CALL(ch_gfa_raw.toList(), ch_assemblies)
    ch_versions = ch_versions.mix(MINIGRAPH_CALL.out.versions)    

    emit:
    gfa = ch_gfa_raw
    fa = ch_fa
    info = ch_info
    bed = MINIGRAPH_CALL.out.bed
    versions = ch_versions

}