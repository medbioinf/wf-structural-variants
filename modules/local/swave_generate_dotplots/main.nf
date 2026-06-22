process SWAVE_GENERATE_DOTPLOTS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    // container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
    //     ? 'https://depot.galaxyproject.org/singularity/mulled-v2-ddb8b80b33a09f54efd9219c18e1d38acfa18bc8:ae02896ffb35dfc564385b2276a1fbf7862567c2-0'
    //     : 'quay.io/biocontainers/mulled-v2-ddb8b80b33a09f54efd9219c18e1d38acfa18bc8:ae02896ffb35dfc564385b2276a1fbf7862567c2-2'}"
    container "quay.io/swave:latest"
    
    input:
    tuple val(meta), path(alt_fasta)
    path(ref_fasta)
    path(gfa_fasta)

    output:
    tuple val(meta), path("*_dotplot_matrices.npz"), emit: dotplot_matrices
    tuple val(meta), path("**.png")                          , emit: pngs, optional: true
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1'), emit: versions, topic: versions

    script:
    def args = task.ext.args ?: ''
    def max_sv_size_param = params.max_sv_size ? "--max_sv_size ${params.max_sv_size}" : ""
    def kmer_size_param = params.kmer_size ? "--kmer_size ${params.kmer_size}" : ""
    def spec_path_param = params.spec_path ? "--spec_path ${params.spec_path}" : ""
    def save_images_flag = params.save_dotplot_images ? "--save_dotplot_images" : ""
    """
    python3 /app/swave/src/generate_dotplots.py \\
        --alt_fasta $alt_fasta \\
        --ref_fasta $ref_fasta \\
        --gfa_fasta $gfa_fasta \\
        --npz_out_prefix "${meta.id}" \\
        --img_out_prefix "${meta.sample}" \\
        $max_sv_size_param \\
        $kmer_size_param \\
        $spec_path_param \\
        $save_images_flag \\
        $args
    """
}