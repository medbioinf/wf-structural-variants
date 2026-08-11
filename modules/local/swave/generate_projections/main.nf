process SWAVE_GENERATE_PROJECTIONS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'docker://jonahkps/panswave:0.1.0'
        : 'docker.io/jonahkps/panswave:0.1.0'}"
    
    input:
    tuple val(meta), path(dotplot_bundles_pkl)

    output:
    tuple val(meta), path("*.projections.pkl.gz"), emit: projections
    tuple val(meta), path("**.png"), emit: pngs, optional: true
    tuple val("${task.process}"), val('swave'), eval("swave-generate-projections --version 2>/dev/null | tail -n1"), emit: versions_swave, topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    swave-generate-projections \\
        --dotplot_bundles_pkl $dotplot_bundles_pkl \\
        --projections_out_prefix "${prefix}" \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip -c > ${prefix}.projections.pkl.gz
    touch stub.projections.png
    """
}