process SWAVE_SPLIT_ALLELES {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    container "quay.io/swave:latest"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.split_*.fa"), emit: splits
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1'), emit: versions, topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def seq_per_split = task.ext.seq_per_split ?: params.seq_per_split
    """
    swave-split-alleles \\
        --fasta ${fasta} \\
        --seq_per_split ${seq_per_split} \\
        --prefix ${prefix} \\
        ${args}
    """
}