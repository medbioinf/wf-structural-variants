process SWAVE_PREDICT {
    tag "${meta.id}"
    label 'process_medium'

    // TODO: handle gpu support

    conda "${moduleDir}/environment.yml"

    container "quay.io/swave:latest"
    
    input:
    tuple val(meta), path(projections_pkl)

    output:
    tuple val(meta), path("*_predictions.pkl.gz"), emit: predictions
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1'), emit: versions, topic: versions

    script:
    def args = task.ext.args ?: ''
    def device_flag = params.use_gpu ? "--device gpu" : "--device cpu"
    def model_param = params.model ? "--model ${params.model}" : ""
    """
    swave-predict \\
        --projections_pkl $projections_pkl \\
        --predictions_out_prefix "${meta.id}" \\
        --cpu_threads ${task.cpus} \\
        $device_flag \\
        $model_param \\
        $args
    """
}