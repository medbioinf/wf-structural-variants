include { MINIGRAPH_CONSTRUCT } from '../../../modules/local/minigraph_construct/main'
include { MINIGRAPH_CALL } from '../../../modules/local/minigraph_call/main'

workflow PANGENOME_GRAPH {

    take:
    ch_reference
    ch_assemblies

    main:
    ch_versions = channel.empty()
    ch_gfa = channel.empty()

    if (!params.gfa) {
        MINIGRAPH_CONSTRUCT(ch_reference, ch_assemblies)
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