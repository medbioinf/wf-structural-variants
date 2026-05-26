process GFATOOLS_GFA2FA {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/gfatools:0.5--h577a1d6_5'
        : 'quay.io/biocontainers/gfatools:0.5--h577a1d6_5'}"

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path("${gfa.baseName}.gfa2fa.fa"), emit: fa
    tuple val("${task.process}"), val('gfatools'), val('0.5'), emit: versions, topic: versions

    script:
    def args = task.ext.args ?: ''
    """
    gfatools \\
        gfa2fa \\
        $args \\
        $gfa \\
        > ${gfa.baseName}.gfa2fa.fa
    """
}