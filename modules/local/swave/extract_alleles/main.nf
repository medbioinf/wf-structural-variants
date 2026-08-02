process SWAVE_EXTRACT_ALLELES {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    // container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
    //     ? 'https://depot.galaxyproject.org/singularity/python:3.12'
    //     : 'quay.io/biocontainers/python:3.12'}"
    container "quay.io/swave:latest"

    input:
    tuple val(meta), path(bed)
    path(gfa_fasta)

    output:
    tuple val(meta), path("${meta.id}_alleles.fa"), emit: fa
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1'), emit: versions, topic: versions

    script:
    def args = task.ext.args ?: ''
    def force_reverse_flag = params.force_reverse ? "--force_reverse" : ""
    def remove_small_flag = params.remove_small ? "--remove_small" : ""
    def min_sv_size_param = params.min_sv_size ? "--min_sv_size ${params.min_sv_size}" : ""
    def spec_snarl_param = params.spec_snarl ? "--spec_snarl ${params.spec_snarl}" : ""
    """
    swave-extract-alleles \\
        --gfa_fasta $gfa_fasta \\
        --bed $bed \\
        --sample_id ${meta.id} \\
        $force_reverse_flag \\
        $remove_small_flag \\
        $min_sv_size_param \\
        $spec_snarl_param \\
        $args \\
        --output ${meta.id}_alleles.fa
    """
}