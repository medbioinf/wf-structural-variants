process SWAVE_EXTRACT_ALLELES {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml" // TODO: update conda path after publishing

    // TODO: publish swave container to quay.io and update the container path
    // container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
    //     ? ''
    //     : ''}"
    container "quay.io/swave:latest"

    input:
    tuple val(meta), path(bed), path(vcf)
    path(gfa_fasta)

    output:
    tuple val(meta), path("*_alleles.fa"), emit: fa
    tuple val("${task.process}"), val('swave'), eval("swave-extract-alleles --version 2>&1"), emit: versions_swave, topic: versions

    script:
    def args = task.ext.args ?: ''
    def bed_param = bed ? "--bed ${bed}" : ""
    def vcf_param = vcf ? "--vcf ${vcf}" : ""
    """
    swave-extract-alleles \\
        --gfa_fasta $gfa_fasta \\
        --sample_id ${meta.id} \\
        ${bed_param} \\
        ${vcf_param} \\
        ${args} \\
        --output_dir .
    """

    stub:
    """
    touch ${meta.id}_alleles.fa
    """
}