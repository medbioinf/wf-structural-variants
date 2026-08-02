include { SAMTOOLS_CAT } from '../../../modules/nf-core/samtools/cat/main'
include { SAMTOOLS_BAM2FQ } from '../../../modules/nf-core/samtools/bam2fq/main'
include { HIFIASM } from '../../../modules/nf-core/hifiasm/main'
include { GFATOOLS_GFA2FA } from '../../../modules/nf-core/gfatools/gfa2fa/main'

workflow LONGREAD_ASSEMBLY {

    take:
    ch_bams // channel: [ meta, bam ]

    main:
    ch_versions = channel.empty()

    ch_bam_files = ch_bams.map { meta, bam_dir ->
        def dir_path = file(bam_dir)
        def bam_files = files("${dir_path}/**.bam").unique()

        if (!bam_files || bam_files.isEmpty()) {
            error "No BAM files found in '${dir_path.toAbsolutePath()}'."
        }

        def meta_new = meta + [ id: meta.id ?: meta.sample, single_end: true ]
        return [ meta_new, bam_files ]
    }

    SAMTOOLS_CAT(ch_bam_files)
    ch_versions = ch_versions.mix(SAMTOOLS_CAT.out.versions_samtools)

    ch_bams_fastq = SAMTOOLS_CAT.out.bam.map { meta, merged_bam ->
        [ meta + [ single_end: true ], merged_bam ]
    }

    SAMTOOLS_BAM2FQ(ch_bams_fastq, false)
    ch_versions = ch_versions.mix(SAMTOOLS_BAM2FQ.out.versions_samtools)

    ch_hifiasm_in = SAMTOOLS_BAM2FQ.out.reads.map { meta, fastq -> 
        [ meta, fastq, [] ] 
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