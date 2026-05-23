include { MINIGRAPH_CONSTRUCT } from '../../../modules/local/minigraph_construct/main'
include { MINIGRAPH_CALL } from '../../../modules/local/minigraph_call/main'

workflow PANGENOME_GRAPH {

    take:
    ch_reference
    ch_assemblies

    main:
    ch_versions = channel.empty()
    ch_gfa = channel.empty()

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

        ch_gfa = MINIGRAPH_CONSTRUCT.out.gfa.map { _meta, gfa -> gfa }
        ch_versions = ch_versions.mix(MINIGRAPH_CONSTRUCT.out.versions)
    } else {
        ch_gfa = channel.fromPath(params.gfa)
    }

    MINIGRAPH_CALL(ch_gfa, ch_assemblies)
    ch_versions = ch_versions.mix(MINIGRAPH_CALL.out.versions)

    emit:
    gfa = ch_gfa
    bed = MINIGRAPH_CALL.out.bed
    versions = ch_versions

}