include { PANSN_FORMAT } from '../../../modules/local/pansn_format/main'
include { MINIGRAPH_PANGENOME } from '../../../modules/local/minigraph/pangenome/main'
include { MINIGRAPH_CALL } from '../../../modules/local/minigraph/call/main'
include { MINIGRAPH_CALL as MINIGRAPH_CALL_REF } from '../../../modules/local/minigraph/call/main'
include { PGGB } from '../../../modules/local/pggb/main'
include { CACTUS_PANGENOME } from '../../../modules/local/cactus/pangenome/main'
include { GFA_SET_REFERENCE } from '../../../modules/local/gfa_set_reference/main'

include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_REF }   from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_PANGENOME } from '../../../modules/nf-core/samtools/faidx/main'
include { VG_DECONSTRUCT } from '../../../modules/nf-core/vg/deconstruct/main'
include { GFATOOLS_GFA2FA } from '../../../modules/nf-core/gfatools/gfa2fa/main'
include { GUNZIP as GUNZIP_PANGENOME_GFA } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PANGENOME_FA } from '../../../modules/nf-core/gunzip/main'


workflow PANGENOME_GRAPH {

    take:
    ch_reference    // channel: [ [id:'reference', sample:'reference', haplotype:0, is_ref:true], fasta ]
    ch_assemblies   // channel: [ meta, fasta ]

    main:
    ch_versions = channel.empty()

    def external_gfa = params.gfa != null
    def ch_reference_adjusted = ch_reference

    // reads gfa file to determine the reference sample name to align the provided reference FASTA's PanSN naming with (in case of differing reference sample name)
    if (external_gfa) {
        def lines = new File(params.gfa).readLines()
        def ref_line = lines.find { line -> line.startsWith('S\t') && line.contains('SR:i:0') }

        def ch_ref_sample = null
        if (ref_line) {
            def match = (ref_line =~ /SN:Z:([^\t]+)/)
            if (match) {
                ch_ref_sample = match[0][1].split('#')[0]
            }
        }

        if (!ch_ref_sample) {
            if (params.graph_construction_tool == "minigraph") {
                error "Could not determine the reference sample name from --gfa (${params.gfa}). Expected a GFA with a reference segment tagged SR:i:0."
            } else {
                log.warn "Could not detect a reference sample name in --gfa (${params.gfa}) (no GFA reference tag, as expected for pggb/cactus GFAs). Ensure your --fasta's derived sample name already matches the reference naming used inside the provided GFA."
            }
        } else {
            ch_reference_adjusted = ch_reference.map { meta, fa ->
                [ meta + [ sample: ch_ref_sample, haplotype: '0' ], fa ]
            }
        }
    }

    ch_all_input_fastas = ch_reference_adjusted.mix(ch_assemblies)
    PANSN_FORMAT(ch_all_input_fastas)
    ch_versions = ch_versions.mix(PANSN_FORMAT.out.versions_gawk)

    PANSN_FORMAT.out.fasta
        .branch { meta, _fa ->
            ref: meta.is_ref == true
            assemblies: meta.is_ref != true
        }
        .set { ch_pansn_split }

    ch_formatted_assemblies = ch_pansn_split.assemblies
    ch_formatted_reference = ch_pansn_split.ref
    ch_ref_fasta_for_swave = ch_formatted_reference

    ch_ref_contig_names = channel.empty()
    if (params.graph_construction_tool == "pggb" || params.graph_construction_tool == "cactus") {
        ch_reference_indexed = ch_formatted_reference.map { meta, fa -> [ meta, fa, [] ] }
        SAMTOOLS_FAIDX_REF(ch_reference_indexed, false)
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_REF.out.versions_samtools)

        ch_ref_contig_names = SAMTOOLS_FAIDX_REF.out.fai
            .map { _meta, fai -> fai }
    }

    if (!params.gfa || (params.gfa && params.graph_construction_tool == "minigraph" && params.minigraph_incremental)) {

        if (params.graph_construction_tool == "minigraph") {

            ch_minigraph_assemblies = ch_formatted_assemblies.map { _meta, fa -> fa }.collect()

            if (params.minigraph_incremental) {
                ch_incremental_ref = channel.fromPath(params.gfa)
                    .map { gfa -> [ [ id: 'pangenome' ], gfa ] }
                MINIGRAPH_PANGENOME(ch_incremental_ref, ch_minigraph_assemblies)
            } else {
                MINIGRAPH_PANGENOME(ch_formatted_reference, ch_minigraph_assemblies)
            }

            ch_gfa_with_meta = MINIGRAPH_PANGENOME.out.gfa
            ch_versions = ch_versions.mix(MINIGRAPH_PANGENOME.out.versions_minigraph)

        } else if (params.graph_construction_tool == "pggb") {

            ch_formatted_reference
                .mix(ch_formatted_assemblies)
                .set { ch_all_fastas }
            
            ch_all_fastas
                .count()
                .set { ch_num_haplotypes }
            
            ch_formatted_reference
                .map { meta, _fa -> meta }
                .combine(ch_num_haplotypes)
                .map { meta, n_count -> meta + [ id: 'pangenome', num_haplotypes: n_count ] }
                .set { ch_graph_meta }
            
            ch_all_fastas
                .map { _meta, fa -> fa }
                .collectFile(name: 'pangenome_input.fa.gz')
                .combine(ch_graph_meta)
                .map { fa, meta -> [ meta, fa, [] ] }
                .set { ch_merged_fasta }
            
            SAMTOOLS_FAIDX_PANGENOME(ch_merged_fasta, false)
            ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_PANGENOME.out.versions_samtools)

            ch_merged_fasta
                .map { meta, fa, _fai -> [ meta, fa ] }
                .join(SAMTOOLS_FAIDX_PANGENOME.out.fai)
                .set { ch_pggb_input }

            PGGB(ch_pggb_input)

            ch_gfa_with_meta = PGGB.out.gfa
            ch_versions = ch_versions.mix(PGGB.out.versions_pggb)

        } else if (params.graph_construction_tool == "cactus") {

            ch_formatted_reference
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

            ch_formatted_reference
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
            ch_versions = ch_versions.mix(CACTUS_PANGENOME.out.versions_cactus)

        } else {
            error "Invalid graph construction tool specified: ${params.graph_construction_tool}. Supported tools are: minigraph, pggb, cactus."
        }

    } else {
        ch_gfa_with_meta = channel.fromPath(params.gfa).map { gfa -> [ [id: gfa.baseName], gfa] }
    }

    GFATOOLS_GFA2FA(ch_gfa_with_meta)
    ch_fa_raw = GFATOOLS_GFA2FA.out.fasta
    ch_versions = ch_versions.mix(GFATOOLS_GFA2FA.out.versions_gfatools)

    ch_fa_raw.branch { _meta, fa ->
        zipped: fa.name.endsWith('.gz')
        unzipped: !fa.name.endsWith('.gz')
    }.set { ch_pangenome_fa_split }

    GUNZIP_PANGENOME_FA(ch_pangenome_fa_split.zipped)
    ch_versions = ch_versions.mix(GUNZIP_PANGENOME_FA.out.versions_gunzip)

    ch_fa = GUNZIP_PANGENOME_FA.out.gunzip
        .mix(ch_pangenome_fa_split.unzipped)
        .map { _meta, fa -> fa }
    ch_gfa_raw = ch_gfa_with_meta.map { _meta, gfa -> gfa }.collect()
    ch_bed = channel.empty()
    ch_ref_bed = channel.empty()
    ch_vcf = channel.empty()

    if (!params.pangenome_only) {
        if (params.graph_construction_tool == "minigraph") {
            MINIGRAPH_CALL_REF(ch_gfa_raw, ch_formatted_reference)
            ch_ref_bed = MINIGRAPH_CALL_REF.out.bed
            ch_versions = ch_versions.mix(MINIGRAPH_CALL_REF.out.versions_minigraph)

            ch_sample_call_gated = ch_formatted_assemblies
                .combine(MINIGRAPH_CALL_REF.out.bed.map { _meta, _bed -> true })
                .map { meta, fa, _gate -> [ meta, fa ] }

            MINIGRAPH_CALL(ch_gfa_raw, ch_sample_call_gated)
            ch_bed = MINIGRAPH_CALL.out.bed
            ch_versions = ch_versions.mix(MINIGRAPH_CALL.out.versions_minigraph)
        } else if (params.graph_construction_tool == "pggb" || params.graph_construction_tool == "cactus") {
            ch_gfa_with_meta
                .branch { _meta, gfa ->
                    zipped: gfa.name.endsWith('.gz')
                    unzipped: !gfa.name.endsWith('.gz')
                }
                .set { ch_pangenome_gfa_split }

            GUNZIP_PANGENOME_GFA(ch_pangenome_gfa_split.zipped)
            ch_versions = ch_versions.mix(GUNZIP_PANGENOME_GFA.out.versions_gunzip)

            ch_decompressed_gfa = GUNZIP_PANGENOME_GFA.out.gunzip
                .mix(ch_pangenome_gfa_split.unzipped)

            GFA_SET_REFERENCE(ch_decompressed_gfa, ch_ref_contig_names.toList())
            ch_versions = ch_versions.mix(GFA_SET_REFERENCE.out.versions_gawk)

            VG_DECONSTRUCT(GFA_SET_REFERENCE.out.gfa, [], [])
            ch_vcf = VG_DECONSTRUCT.out.vcf
            ch_versions = ch_versions.mix(VG_DECONSTRUCT.out.versions_vg)
        }
    }

    emit:
    pangenome_fa = ch_fa
    ref_fasta_zipped = ch_ref_fasta_for_swave
    bed = ch_bed
    ref_bed = ch_ref_bed
    vcf = ch_vcf
    versions = ch_versions
}