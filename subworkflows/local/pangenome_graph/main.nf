include { MINIGRAPH } from '../../../modules/local/minigraph/main'

workflow PANGENOME_GRAPH {

    take:
    ch_reference
    ch_assemblies

    main:
    ch_versions = channel.empty()

    MINIGRAPH(ch_reference, ch_assemblies)

    ch_versions = ch_versions.mix(MINIGRAPH.out.versions)

    emit:
    gfa = MINIGRAPH.out.gfa
    versions = ch_versions

}