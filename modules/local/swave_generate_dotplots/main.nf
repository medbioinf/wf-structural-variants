process SWAVE_GENERATE_DOTPLOTS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/matplotlib:3.5.1'
        : 'quay.io/biocontainers/matplotlib:3.5.1'}"
    
    input:
    tuple val(meta), path(alleles_fasta)
    path(gfa_fasta)

    output:
    tuple val(meta), path("${meta.id}_matrices.npz"), emit: matrices
    tuple val(meta), path("*.png"), emit: pngs, optional: true
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1'), emit: versions, topic: versions

    script:
    def args = task.ext.args ?: ''
    def save_images_flag = params.save_dotplot_images ? "--save_dotplot_images" : ""
    def kmer_size_param = params.kmer_size ? "--kmer_size ${params.kmer_size}" : ""
    def max_dim_param = params.max_dotplot_dim ? "--max_dotplot_dim ${params.max_dotplot_dim}" : ""
    """
    generate_dotplots.py \\
        --alleles_fasta $alleles_fasta \\
        --gfa_fasta $gfa_fasta \\
        --out_prefix "${meta.id}" \\
        $save_images_flag \\
        $kmer_size_param \\
        $max_dim_param \\
        $args
    """
}