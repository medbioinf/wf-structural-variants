process SWAVE_EXTRACT_VARIANT_SEQUENCES {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'docker://jonahkps/panswave:0.1.0'
        : 'docker.io/jonahkps/panswave:0.1.0'}"

    input:
    tuple val(meta), path(vcf)
    path(gfa_fasta)

    output:
    tuple val(meta), path("*.variant_sequences.fa"), emit: fa
    tuple val("${task.process}"), val('swave'), eval("swave-extract-variant-sequences --version 2>/dev/null | tail -n1"), emit: versions_swave, topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    swave-extract-variant-sequences \\
        --vcf $vcf \\
        --gfa_fasta $gfa_fasta \\
        --output_fasta ${prefix}.variant_sequences.fa \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.variant_sequences.fa
    """
}