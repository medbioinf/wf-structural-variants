process SWAVE_SPLIT_ALLELES {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml" // TODO: update conda path after publishing

    // TODO: publish swave container to quay.io and update the container path
    // container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
    //     ? ''
    //     : ''}"
    container "quay.io/swave:latest"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.split_*.fa"), emit: splits
    tuple val("${task.process}"), val('swave'), eval("swave-split-alleles --version 2>/dev/null | tail -n1"), emit: versions_swave, topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def seq_per_split = task.ext.seq_per_split ?: 2500
    """
    swave-split-alleles \\
        --fasta ${fasta} \\
        --seq_per_split ${seq_per_split} \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.split_01.fa
    """
}