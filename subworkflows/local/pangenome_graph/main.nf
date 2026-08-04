include { PANSN_FORMAT }        from '../../../modules/local/pansn_format/main'
include { MINIGRAPH_CONSTRUCT } from '../../../modules/local/minigraph/construct/main'
include { MINIGRAPH_CALL }      from '../../../modules/local/minigraph/call/main'
include { PGGB }                from '../../../modules/local/pggb/main'
include { CACTUS_PANGENOME }    from '../../../modules/local/cactus/pangenome/main'
include { GFA_SET_REFERENCE }   from '../../../modules/local/gfa_set_reference/main'

include { SAMTOOLS_FAIDX }      from '../../../modules/nf-core/samtools/faidx/main'
include { VG_DECONSTRUCT }      from '../../../modules/nf-core/vg/deconstruct/main'
include { GFATOOLS_GFA2FA }     from '../../../modules/nf-core/gfatools/gfa2fa/main'
include { GUNZIP as GUNZIP_FA }  from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFA } from '../../../modules/nf-core/gunzip/main'

workflow PANGENOME_GRAPH {

    take:
    ch_reference   // channel: [ [id:'reference', sample:'reference', haplotype:0, is_ref:true], fasta ]
    ch_assemblies  // channel: [ meta, fasta ]

    main:
    ch_versions = channel.empty()

    ch_all_input_fastas = ch_reference.mix(ch_assemblies)
    
    PANSN_FORMAT(ch_all_input_fastas)
    ch_versions = ch_versions.mix(PANSN_FORMAT.out.versions_gawk)

    PANSN_FORMAT.out.fasta
        .branch { meta, _fa ->
            ref: meta.is_ref == true
            assemblies: meta.is_ref != true
        }
        .set { ch_pansn_split }

    ch_formatted_ref = ch_pansn_split.ref
    ch_formatted_assemblies = ch_pansn_split.assemblies

    if (!params.gfa) {

        if (params.graph_construction_tool == "minigraph") {

            ch_minigraph_assemblies = ch_formatted_assemblies.map { _meta, fa -> fa }.collect()
            MINIGRAPH_CONSTRUCT(ch_formatted_ref, ch_minigraph_assemblies)

            ch_gfa_with_meta = MINIGRAPH_CONSTRUCT.out.gfa
            ch_info = MINIGRAPH_CONSTRUCT.out.info  // TODO: remove
            ch_versions = ch_versions.mix(MINIGRAPH_CONSTRUCT.out.versions)

        } else if (params.graph_construction_tool == "pggb") {

            ch_formatted_ref
                .mix(ch_formatted_assemblies)
                .set { ch_all_fastas }
            
            ch_all_fastas
                .count()
                .set { ch_num_haplotypes }
            
            ch_formatted_ref
                .map { meta, _fa -> meta }
                .combine(ch_num_haplotypes)
                .map { meta, n_count -> meta + [ id: 'pangenome', num_haplotypes: n_count ] }
                .set { ch_graph_meta }
            
            ch_all_fastas
                .map { _meta, fa -> fa }
                .collectFile(name: 'pangenome_input.fasta', newLine: true)
                .combine(ch_graph_meta)
                .map { fa, meta -> [ meta, fa, [] ] }
                .set { ch_merged_fasta }
            
            SAMTOOLS_FAIDX(ch_merged_fasta, false)
            ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions_samtools)

            ch_merged_fasta
                .map { meta, fa, _fai -> [ meta, fa ] }
                .join(SAMTOOLS_FAIDX.out.fai)
                .set { ch_pggb_input }

            PGGB(ch_pggb_input)

            ch_gfa_with_meta = PGGB.out.gfa
            ch_info = channel.empty()
            ch_versions = ch_versions.mix(PGGB.out.versions_pggb)

        } else if (params.graph_construction_tool == "cactus") {

            ch_formatted_ref
                .mix(ch_formatted_assemblies)
                .set { ch_all_fastas_cactus }

            ch_seqfile = ch_all_fastas_cactus
                .map { meta, fa ->
                    def sample_safe = meta.sample.replaceAll(/\./, '_')
                    def seq_name = meta.is_ref ? sample_safe : "${sample_safe}.${meta.haplotype}"
                    "${seq_name}\t${fa.name}"
                }
                .collectFile(name: 'seqfile.txt', newLine: true, sort: false)

            ch_num_haplotypes_cactus = ch_all_fastas_cactus.count()

            ch_formatted_ref
                .map { meta, _fa -> meta }
                .combine(ch_num_haplotypes_cactus)
                .map { meta, n_count ->
                    meta + [ id: 'pangenome', sample: meta.sample.replaceAll(/\./, '_'), num_haplotypes: n_count ]
                }
                .combine(ch_seqfile)
                .set { ch_cactus_meta_seqfile }

            ch_cactus_fastas = ch_all_fastas_cactus
                .map { _meta, fa -> fa }
                .collect()

            CACTUS_PANGENOME(ch_cactus_meta_seqfile, ch_cactus_fastas)

            ch_gfa_with_meta = CACTUS_PANGENOME.out.gfa
            ch_info = channel.empty()
            ch_versions = ch_versions.mix(CACTUS_PANGENOME.out.versions_cactus)
            
        } else {
            error "Invalid graph construction tool specified: ${params.graph_construction_tool}. Supported tools are: minigraph, pggb, cactus."
        }

    } else {
        ch_gfa_with_meta = channel.fromPath(params.gfa).map { gfa -> [ [id: gfa.baseName], gfa] }
    }

    if (params.gfa2fa_fa) { // TODO: potentially unnecessary, remove?
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

    GUNZIP_FA(ch_fa_split.zipped)
    ch_versions = ch_versions.mix(GUNZIP_FA.out.versions_gunzip)

    ch_fa = GUNZIP_FA.out.gunzip
        .mix(ch_fa_split.unzipped)
        .map { _meta, fa -> fa }
    
    ch_gfa_raw = ch_gfa_with_meta.map { _meta, gfa -> gfa }.collect()
    ch_bed = channel.empty()
    ch_vcf = channel.empty()

    if (params.graph_construction_tool == "minigraph") {
        ch_all_call_inputs = ch_formatted_ref.mix(ch_formatted_assemblies)

        MINIGRAPH_CALL(ch_gfa_raw, ch_all_call_inputs)
        ch_bed = MINIGRAPH_CALL.out.bed
        ch_versions = ch_versions.mix(MINIGRAPH_CALL.out.versions)
    } else if (params.graph_construction_tool == "pggb" || params.graph_construction_tool == "cactus") {
        ch_gfa_with_meta
            .branch { _meta, gfa ->
                zipped: gfa.name.endsWith('.gz')
                unzipped: !gfa.name.endsWith('.gz')
            }
            .set { ch_gfa_deconstruct_split }

        GUNZIP_GFA(ch_gfa_deconstruct_split.zipped)
        ch_versions = ch_versions.mix(GUNZIP_GFA.out.versions_gunzip)

        ch_decompressed_gfa = GUNZIP_GFA.out.gunzip
            .mix(ch_gfa_deconstruct_split.unzipped)
        
        GFA_SET_REFERENCE(ch_decompressed_gfa)
        ch_versions = ch_versions.mix(GFA_SET_REFERENCE.out.versions_gawk)

        VG_DECONSTRUCT(GFA_SET_REFERENCE.out.gfa, [], [])
        ch_vcf = VG_DECONSTRUCT.out.vcf
        ch_versions = ch_versions.mix(VG_DECONSTRUCT.out.versions_vg)
    }

    emit:
    gfa = ch_gfa_raw
    pangenome_fa = ch_fa
    ref_fasta_pansn = ch_formatted_ref.map { _meta, fa -> fa }
    info = ch_info
    bed = ch_bed
    vcf = ch_vcf
    versions = ch_versions
}