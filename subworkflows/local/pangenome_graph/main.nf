include { MINIGRAPH_CONSTRUCT } from '../../../modules/local/minigraph/construct/main'
include { MINIGRAPH_CALL } from '../../../modules/local/minigraph/call/main'
include { GFATOOLS_GFA2FA } from '../../../modules/nf-core/gfatools/gfa2fa/main'
include { GUNZIP } from '../../../modules/nf-core/gunzip/main'

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
            .branch { meta, _fa ->
                samples: meta.is_ref != true
                ref:    meta.is_ref == true
            }
            .set { ch_split_for_construct }
        
        ch_split_for_construct.samples
            .multiMap { meta, fa ->
                metas: meta
                fastas: fa
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

    if (params.gfa2fa_fa) {
        ch_fa_raw = channel.fromPath(params.gfa2fa_fa).map { fa -> [ [id: fa.baseName], fa ] }
    } else {
        GFATOOLS_GFA2FA(ch_gfa_with_meta)
        ch_fa_raw = GFATOOLS_GFA2FA.out.fasta
        ch_versions = ch_versions.mix(GFATOOLS_GFA2FA.out.versions_gfatools)
    }

    ch_fa_raw.branch { _meta, fa ->
        zipped: fa.name.endsWith('.gz')
        unzipped: !fa.name.endsWith('.gz')
    }.set { ch_fa_split }

    GUNZIP(ch_fa_split.zipped)
    ch_versions = ch_versions.mix(GUNZIP.out.versions_gunzip)

    ch_fa = GUNZIP.out.gunzip
        .mix(ch_fa_split.unzipped)
        .map { _meta, fa -> fa }

    ch_gfa_raw = ch_gfa_with_meta.map { _meta, gfa -> gfa }.collect()
    MINIGRAPH_CALL(ch_gfa_raw, ch_assemblies)
    ch_versions = ch_versions.mix(MINIGRAPH_CALL.out.versions)    

    emit:
    gfa = ch_gfa_raw
    fa = ch_fa
    info = ch_info
    bed = MINIGRAPH_CALL.out.bed
    versions = ch_versions
}