process SWAVE_PREDICT {
    tag "${meta.id}"
    label 'process_medium'

    // TODO: handle gpu support

    conda "${moduleDir}/environment.yml" // TODO: update conda path after publishing

    // TODO: publish swave container to quay.io and update the container path
    // container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
    //     ? ''
    //     : ''}"
    container "quay.io/swave:latest"
    
    input:
    tuple val(meta), path(projections_pkl)

    output:
    tuple val(meta), path("*.predictions.pkl.gz"), emit: predictions
    tuple val("${task.process}"), val('swave'), eval("swave-predict --version 2>/dev/null | tail -n1"), emit: versions_swave, topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    swave-predict \\
        --projections_pkl $projections_pkl \\
        --predictions_out_prefix "${prefix}" \\
        --cpu_threads ${task.cpus} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip -c > ${prefix}.predictions.pkl.gz
    """
}