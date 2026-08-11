include { SWAVE_EXTRACT_ALLELES } from '../../../modules/local/swave/extract_alleles/main'
include { SWAVE_SPLIT_ALLELES } from '../../../modules/local/swave/split_alleles/main'
include { SWAVE_GENERATE_DOTPLOTS } from '../../../modules/local/swave/generate_dotplots/main'
include { SWAVE_GENERATE_PROJECTIONS } from '../../../modules/local/swave/generate_projections/main'

include { BCFTOOLS_QUERY as BCFTOOLS_QUERY_LIST_SAMPLES } from '../../../modules/nf-core/bcftools/query/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_SAMPLES } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_REF } from '../../../modules/nf-core/bcftools/view/main'

workflow SWAVE_PREPROCESSING {

    take:
    ch_bed          // channel: [ meta, bed ]
    ch_vcf          // channel: [ meta, vcf ]
    ch_gfa_fasta    // channel: [ fasta ]
    ch_ref_fasta    // channel: [ fasta ]

    main:
    ch_versions = channel.empty()

    if (params.graph_construction_tool == "minigraph") {

        ch_extract_input = ch_bed.map { meta, bed -> [ meta, meta.is_ref, bed, [] ] }

        ch_ref_bed = ch_bed
            .filter { meta, _bed -> meta.is_ref }
            .map { _meta, bed -> bed }

        SWAVE_EXTRACT_ALLELES(ch_extract_input, ch_ref_bed.toList(), ch_gfa_fasta.toList())
        ch_versions = ch_versions.mix(SWAVE_EXTRACT_ALLELES.out.versions_swave)
        ch_extracted_alleles = SWAVE_EXTRACT_ALLELES.out.fa
        ch_extracted_equal_paths = SWAVE_EXTRACT_ALLELES.out.equal_paths

    } else {

        ch_vcf_indexed = ch_vcf.map { meta, vcf -> [ meta, vcf, [] ] }

        BCFTOOLS_VIEW_REF(ch_vcf_indexed, [], [], [])
        ch_versions = ch_versions.mix(BCFTOOLS_VIEW_REF.out.versions_bcftools)

        ch_ref_vcf = BCFTOOLS_VIEW_REF.out.vcf
            .map { meta, vcf -> [ meta + [ id: "${meta.sample}", is_ref: true ], vcf ] }
        
        BCFTOOLS_QUERY_LIST_SAMPLES(ch_vcf_indexed, [], [], [])
        ch_versions = ch_versions.mix(BCFTOOLS_QUERY_LIST_SAMPLES.out.versions_bcftools)

        ch_sample_names = BCFTOOLS_QUERY_LIST_SAMPLES.out.output
            .flatMap { _meta, txt -> txt.readLines() }

        ch_vcf_val = ch_vcf.map { _meta, vcf -> vcf }

        ch_vcf_per_sample = ch_sample_names
            .combine(ch_vcf_val)
            .map { sample_name, vcf -> [ [ id: sample_name, sample: sample_name, is_ref: false ], vcf, [] ] }

        BCFTOOLS_VIEW_SAMPLES(ch_vcf_per_sample, [], [], [])
        ch_versions = ch_versions.mix(BCFTOOLS_VIEW_SAMPLES.out.versions_bcftools)

        ch_vcf_for_extraction = ch_ref_vcf.mix(BCFTOOLS_VIEW_SAMPLES.out.vcf)

        ch_extract_input = ch_vcf_for_extraction.map { meta, vcf -> [ meta, false, [], vcf ] }

        SWAVE_EXTRACT_ALLELES(ch_extract_input, channel.empty().toList(), ch_gfa_fasta.toList())
        ch_versions = ch_versions.mix(SWAVE_EXTRACT_ALLELES.out.versions_swave)
        ch_extracted_alleles = SWAVE_EXTRACT_ALLELES.out.fa
        ch_extracted_equal_paths = SWAVE_EXTRACT_ALLELES.out.equal_paths
    }

    ch_extracted_alleles = ch_extracted_alleles
        .flatMap { meta, fa_files ->
            def files = fa_files instanceof List ? fa_files : [fa_files]
            files.collect { fa ->
                def label = fa.baseName - '_alleles'
                def new_meta = meta + [ id: label ]
                return [ new_meta, fa ]
            }
        }

    ch_extracted_equal_paths = ch_extracted_equal_paths
        .flatMap { meta, txt_files ->
            def files = txt_files instanceof List ? txt_files : [txt_files]
            files.collect { txt ->
                def label = txt.baseName - '_equal_paths'
                def new_meta = meta + [ id: label ]
                return [ new_meta, txt ]
            }
        }

    ch_sample_alleles = ch_extracted_alleles
        .filter { meta, _fa -> meta.is_ref != true }
    SWAVE_SPLIT_ALLELES(ch_sample_alleles)
    ch_versions = ch_versions.mix(SWAVE_SPLIT_ALLELES.out.versions_swave)

    ch_equal_paths_combined = ch_extracted_equal_paths
        .filter { meta, _txt -> meta.is_ref != true }
        .collectFile(name: 'equal_paths_combined.txt', sort: false) { meta, txt_file ->
            txt_file.readLines().collect { line -> "${meta.id}\t${line}\n" }.join('')
        }

    ch_dotplot_inputs = SWAVE_SPLIT_ALLELES.out.splits
        .transpose()
        .map { meta, fa ->
            def new_meta = meta.clone()
            new_meta.sample = meta.id
            new_meta.id = fa.baseName
            return [ new_meta, fa ]
        }

    SWAVE_GENERATE_DOTPLOTS(ch_dotplot_inputs, ch_ref_fasta.toList(), ch_gfa_fasta.toList())
    ch_versions = ch_versions.mix(SWAVE_GENERATE_DOTPLOTS.out.versions_swave)
    ch_dotplots = SWAVE_GENERATE_DOTPLOTS.out.dotplots

    SWAVE_GENERATE_PROJECTIONS(SWAVE_GENERATE_DOTPLOTS.out.dotplots)
    ch_versions = ch_versions.mix(SWAVE_GENERATE_PROJECTIONS.out.versions_swave)
    ch_projections = SWAVE_GENERATE_PROJECTIONS.out.projections

    emit:
    equal_paths = ch_equal_paths_combined
    dotplots = ch_dotplots
    projections = ch_projections
    versions = ch_versions
}