process SWAVE_GENERATE_DOTPLOTS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/matplotlib:3.5.1'
        : 'quay.io/biocontainers/matplotlib:3.5.1'}"
    
    input:
    tuple val(meta), val(ref_seq), val(alt_seq)

    output:
    tuple val(meta), path("${meta.id}_${meta.snarl}_matrices.npz"), emit: matrices
    tuple val(meta), path("*.png"), emit: pngs, optional: true
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1'), emit: versions, topic: versions

    script:
    def args = task.ext.args ?: ''
    def save_images_flag = params.save_dotplot_images ? "--save_dotplot_images" : ""
    def kmer_size_param = params.kmer_size ? "--kmer_size ${params.kmer_size}" : ""
    def max_dim_param = params.max_dotplot_dim ? "--max_dotplot_dim ${params.max_dotplot_dim}" : ""
    def output_prefix = "${meta.id}_${meta.snarl}"
    """
    generate_dotplots.py \\
        --ref_seq "${ref_seq}" \\
        --alt_seq "${alt_seq}" \\
        --out_prefix "${output_prefix}" \\
        $save_images_flag \\
        $kmer_size_param \\
        $max_dim_param \\
        $args
    """
}