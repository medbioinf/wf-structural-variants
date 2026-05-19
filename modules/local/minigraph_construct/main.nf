process MINIGRAPH_CONSTRUCT {
    tag "${meta.id}"
    label 'process_high_single_task'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/minigraph:0.20--he4a0461_2'
        : 'quay.io/biocontainers/minigraph:0.20--he4a0461_2'}"

    input:
    tuple val(meta), path(reference)
    path assemblies

    output:
    tuple val(meta), path("${meta.id}_pangenome.gfa"), emit: gfa
    tuple val("${task.process}"), val('minigraph'), eval('minigraph --version 2>&1'), emit: versions, topic: versions

    script:
    def args = task.ext.args ?: ''
    """
    minigraph \\
        -cxggs \\
        -t $task.cpus \\
        $args \\
        $reference \\
        $assemblies \\
        > ${meta.id}_pangenome.gfa
    """
}