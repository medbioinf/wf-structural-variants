process SWAVE_CALL_VARIANTS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "quay.io/swave:latest"
    
    input:
    tuple val(meta), path(predictions_pkl), path(projections_pkl), path(dotplots_pkl)

    output:
    tuple val(meta), path("*_variants.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1'), emit: versions, topic: versions

    script:
    def args = task.ext.args ?: ''
    def min_sv_size  = params.min_sv_size ? "--min_sv_size ${params.min_sv_size}" : ""
    def max_sv_size  = params.max_sv_size ? "--max_sv_size ${params.max_sv_size}" : ""
    def max_sv_comps = params.max_sv_comps ? "--max_sv_comps ${params.max_sv_comps}" : ""
    def dup_to_ins   = params.dup_to_ins ? "--dup_to_ins" : ""
    """
    python3 /app/swave/src/call_variants.py \\
        --predictions_pkl $predictions_pkl \\
        --projections_pkl $projections_pkl \\
        --dotplots_pkl $dotplots_pkl \\
        --output_tsv "${meta.id}_variants.tsv" \\
        --sample_id "${meta.sample}" \\
        $min_sv_size \\
        $max_sv_size \\
        $max_sv_comps \\
        $dup_to_ins \\
        $args
    """
}
