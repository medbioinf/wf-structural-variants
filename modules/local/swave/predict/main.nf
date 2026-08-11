process SWAVE_PREDICT {
    tag "${meta.id}"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'docker://jonahkps/panswave:0.1.0'
        : 'docker.io/jonahkps/panswave:0.1.0'}"
    
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