process SWAVE_GENERATE_DOTPLOTS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml" // TODO: update conda path after publishing

    // TODO: publish swave container to quay.io and update the container path
    // container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
    //     ? ''
    //     : ''}"
    container "quay.io/swave:latest"
    
    input:
    tuple val(meta), path(alt_fasta)
    path(ref_fasta)
    path(gfa_fasta)

    output:
    tuple val(meta), path("*_dotplots.pkl.gz"), emit: dotplots
    tuple val(meta), path("**.png"), emit: pngs, optional: true
    tuple val("${task.process}"), val('swave'), eval("swave-generate-dotplots --version 2>/dev/null | tail -n1"), emit: versions_swave, topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${meta.sample}
    swave-generate-dotplots \\
        --alt_fasta $alt_fasta \\
        --ref_fasta $ref_fasta \\
        --gfa_fasta $gfa_fasta \\
        --pkl_out_prefix "${prefix}" \\
        --img_out_prefix "${meta.sample}" \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${meta.sample}
    echo "" | gzip -c > ${prefix}_dotplots.pkl.gz
    touch ${meta.sample}/stub_dotplots.png
    """
}