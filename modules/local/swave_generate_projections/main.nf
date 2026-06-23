process SWAVE_GENERATE_PROJECTIONS  {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    container "quay.io/swave:latest"
    
    input:
    tuple val(meta), path(dotplot_bundles_pkl)

    output:
    tuple val(meta), path("*_projections.pkl.gz"), emit: projections
    tuple val(meta), path("**.png"), emit: pngs, optional: true
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1'), emit: versions, topic: versions

    script:
    def args = task.ext.args ?: ''
    def kmer_size_param = params.kmer_size ? "--kmer_size ${params.kmer_size}" : ""
    def save_images_flag = params.save_projections_images ? "--save_projections_images" : ""
    """
    python3 /app/swave/src/generate_projections.py \\
        --dotplot_bundles_pkl $dotplot_bundles_pkl \\
        --projections_out_prefix "${meta.id}" \\
        $kmer_size_param \\
        $save_images_flag \\
        $args
    """
}