include { SAMTOOLS_FASTQ } from '../../../modules/nf-core/samtools/fastq/main'
include { CAT_FASTQ } from '../../../modules/nf-core/cat/fastq/main'
include { HIFIASM } from '../../../modules/nf-core/hifiasm/main'
include { GFATOOLS_GFA2FA } from '../../../modules/nf-core/gfatools/gfa2fa/main'

workflow LONGREAD_ASSEMBLY {

    take:
    ch_bams // channel: [ meta, bam ]

    main:
    ch_versions = channel.empty()

    ch_individual_bams = ch_bams.flatMap { meta, bam_dir ->
        def dir_path = file(bam_dir)
        def bam_files = files("${dir_path}/**.bam").unique()

        if (!bam_files || bam_files.isEmpty()) {
            error "No BAM files found in '${dir_path.toAbsolutePath()}'."
        }

        def sample_id = meta.id ?: meta.sample

        return bam_files.collect { bam ->
            def bam_name = bam.name.replaceAll(/\.bam$/, '')
            def meta_new = meta + [ sample_id: sample_id, id: "${sample_id}_${bam_name}", single_end: true ]
            [ meta_new, bam ]
        }
    }

    SAMTOOLS_FASTQ(ch_individual_bams, false)
    ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions_samtools)

    ch_fastqs_to_cat = SAMTOOLS_FASTQ.out.other
        .map { meta, fastq -> 
            def meta_clean = meta + [ id: meta.sample_id, single_end: true ]
            [ meta_clean, fastq ] 
        }
        .groupTuple()
        .map { meta, fastqs -> 
            [ meta, fastqs.flatten() ] 
        }

    CAT_FASTQ(ch_fastqs_to_cat)
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions_cat)

    ch_hifiasm_in = CAT_FASTQ.out.reads.map { meta, merged_fastq -> 
        [ meta, merged_fastq, [] ] 
    }

    HIFIASM ( 
        ch_hifiasm_in,
        [ [id:'empty_trio'], [], [] ],
        [ [id:'empty_hic'], [], [] ],
        [ [id:'empty_bin'], [] ]
    )
    ch_versions = ch_versions.mix(HIFIASM.out.versions_hifiasm)

    ch_gfa_hap1 = HIFIASM.out.hap1_contigs.map { meta, gfa -> 
        [ meta + [ haplotype: 1, id: "${meta.sample}_hap1", is_ref: false ], gfa ] 
    }
    ch_gfa_hap2 = HIFIASM.out.hap2_contigs.map { meta, gfa -> 
        [ meta + [ haplotype: 2, id: "${meta.sample}_hap2", is_ref: false ], gfa ] 
    }
    
    ch_gfas = ch_gfa_hap1.mix(ch_gfa_hap2)

    GFATOOLS_GFA2FA (ch_gfas)
    ch_versions = ch_versions.mix(GFATOOLS_GFA2FA.out.versions_gfatools)

    emit:
    assemblies = GFATOOLS_GFA2FA.out.fasta
    versions = ch_versions
}