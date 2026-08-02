/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_pangenomesv_pipeline'

include { LONGREAD_ASSEMBLY } from '../subworkflows/local/longread_assembly/main'
include { PANGENOME_GRAPH } from '../subworkflows/local/pangenome_graph/main'
include { SWAVE_PREPROCESSING } from '../subworkflows/local/swave_preprocessing/main'
include { SWAVE_GENOTYPING } from '../subworkflows/local/swave_genotyping/main'
include { SWAVE_ANNOTATION } from '../subworkflows/local/swave_annotation/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PANGENOMESV {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:
    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()


    //
    // Samplesheet processing & optionally run Subworkflow LONGREAD_ASSEMBLY
    //
    ch_samplesheet
        .branch { meta, fasta, bam_dir ->
            fasta: fasta
                return [ meta, file(fasta) ]
            bam: bam_dir
                return [ meta, file(bam_dir) ]
        }
        .set { ch_input_branched }
    
    ch_existing_assemblies = ch_input_branched.fasta
        .map { meta, fasta ->
            def new_meta = meta + [ id: "${meta.sample}_hap${meta.haplotype}", is_ref: false ]
            return [ new_meta, fasta ]
        }
    
    LONGREAD_ASSEMBLY(ch_input_branched.bam)
    ch_assembled_fastas = LONGREAD_ASSEMBLY.out.assemblies

    ch_assemblies = ch_existing_assemblies.mix( ch_assembled_fastas )


    //
    // SUBWORKFLOW: Run Minigraph Pangenome Graph Construction & Snarl Calling
    //
    ch_reference = channel.fromPath(params.fasta, checkIfExists: true).map { fasta -> [ [ id:'reference' ], fasta ] }    
    
    ch_ref_as_assembly = ch_reference.map { meta, fasta -> 
        [ meta + [ sample: 'reference', haplotype: 0, id: 'reference', is_ref: true ], fasta ] 
    }

    ch_assemblies = ch_assemblies.mix(ch_ref_as_assembly)

    PANGENOME_GRAPH(ch_reference, ch_assemblies)

    ch_pangenome_fa = PANGENOME_GRAPH.out.fa
    ch_bed_files = PANGENOME_GRAPH.out.bed
        .filter { meta,_bed -> !meta.is_ref }


    //
    // SUBWORKFLOW: Run SWAVE Preprocessing
    //
    SWAVE_PREPROCESSING(ch_bed_files, ch_pangenome_fa, ch_reference.map{ _meta, fa -> fa })

    ch_dotplots = SWAVE_PREPROCESSING.out.dotplots
    ch_projections = SWAVE_PREPROCESSING.out.projections


    //
    // SUBWORKFLOW: Run SWAVE Genotyping
    //
    SWAVE_GENOTYPING(ch_dotplots, ch_projections, ch_reference.map{ _meta, fa -> fa })


    //
    // SUBWORKFLOW: Run SWAVE Annotation
    //
    ch_vcf_for_annotation = SWAVE_GENOTYPING.out.vcf_split
        .filter { _meta, vcf ->
            vcf.exists() && vcf.readLines().any { line -> !line.startsWith('#') && line.trim() }    // check if VCF has any non-header lines
        }
        
    SWAVE_ANNOTATION ( ch_vcf_for_annotation )


    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'pangenomesv_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'pangenomesv'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
