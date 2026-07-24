process MINIGRAPH_CALL {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/minigraph:0.20--he4a0461_2'
        : 'quay.io/biocontainers/minigraph:0.20--he4a0461_2'}"

    input:
    path(gfa)
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${assembly}.bed"), emit: bed
    tuple val("${task.process}"), val('minigraph'), eval('minigraph --version 2>&1'), emit: versions, topic: versions

    script:
    def args = task.ext.args ?: ''
    """
    minigraph \\
        -t $task.cpus \\
        -xasm \\
        --call \\
        $args \\
        $gfa \\
        $assembly \\
        > ${assembly}.bed
    """
}