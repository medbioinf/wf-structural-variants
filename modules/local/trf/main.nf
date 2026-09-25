process TRF {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/trf:4.09.1--h7b50bb2_7'
        : 'quay.io/biocontainers/trf:4.09.1--h7b50bb2_7'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.dat"), emit: dat
    tuple val("${task.process}"), val('trf'), eval("trf --version 2>&1 | head -n1"), emit: versions_trf, topic: versions

    script:
    def args = task.ext.args ?: '2 7 7 80 10 50 500'
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    trf ${fasta} ${args} -d -h > /dev/null || true
    mv ${fasta}.*.dat ${prefix}.dat
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    touch ${prefix}.dat
    """
}