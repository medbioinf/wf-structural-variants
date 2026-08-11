process SWAVE_CALL_VARIANTS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'docker://jonahkps/panswave:0.1.0'
        : 'docker.io/jonahkps/panswave:0.1.0'}"
    
    input:
    tuple val(meta), path(predictions_pkl), path(projections_pkl), path(dotplots_pkl)

    output:
    tuple val(meta), path("*.variants.tsv"), emit: tsv
    tuple val("${task.process}"), val('swave'), eval("swave-call-variants --version 2>/dev/null | tail -n1"), emit: versions_swave, topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    swave-call-variants \\
        --predictions_pkl $predictions_pkl \\
        --projections_pkl $projections_pkl \\
        --dotplots_pkl $dotplots_pkl \\
        --output_tsv "${prefix}.variants.tsv" \\
        --sample_id "${meta.sample}" \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.variants.tsv
    """
}