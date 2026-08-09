include { SWAVE_EXTRACT_ALLELES } from '../../../modules/local/swave/extract_alleles/main'
include { SWAVE_SPLIT_ALLELES } from '../../../modules/local/swave/split_alleles/main'
include { SWAVE_GENERATE_DOTPLOTS } from '../../../modules/local/swave/generate_dotplots/main'
include { SWAVE_GENERATE_PROJECTIONS } from '../../../modules/local/swave/generate_projections/main'

include { BCFTOOLS_QUERY as BCFTOOLS_QUERY_LIST_SAMPLES } from '../../../modules/nf-core/bcftools/query/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_SAMPLES } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_REF } from '../../../modules/nf-core/bcftools/view/main'

workflow SWAVE_PREPROCESSING {

    take:
    ch_bed
    ch_vcf
    ch_gfa_fasta
    ch_ref_fasta

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

    ch_sample_alleles = ch_extracted_alleles
        .filter { meta, _fa -> meta.is_ref != true }
    
    SWAVE_SPLIT_ALLELES(ch_sample_alleles)
    ch_versions = ch_versions.mix(SWAVE_SPLIT_ALLELES.out.versions_swave)

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
    ch_dotplot_pngs = SWAVE_GENERATE_DOTPLOTS.out.pngs

    SWAVE_GENERATE_PROJECTIONS(SWAVE_GENERATE_DOTPLOTS.out.dotplots)
    ch_versions = ch_versions.mix(SWAVE_GENERATE_PROJECTIONS.out.versions_swave)
    ch_projections = SWAVE_GENERATE_PROJECTIONS.out.projections
    ch_projections_pngs = SWAVE_GENERATE_PROJECTIONS.out.pngs

    emit:
    alleles_fasta = SWAVE_SPLIT_ALLELES.out.splits
    dotplots = ch_dotplots
    dotplot_pngs = ch_dotplot_pngs
    projections = ch_projections
    projections_pngs = ch_projections_pngs
    versions = ch_versions
}