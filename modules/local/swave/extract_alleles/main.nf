process SWAVE_EXTRACT_ALLELES {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'docker://jonahkps/panswave:0.1.0'
        : 'docker.io/jonahkps/panswave:0.1.0'}"

    input:
    tuple val(meta), val(is_ref), path(bed), path(vcf)
    path(ref_bed, stageAs: "ref_bed/*")
    path(gfa_fasta)

    output:
    tuple val(meta), path("*_alleles.fa"), emit: fa
    tuple val(meta), path("*_equal_paths.txt"), emit: equal_paths
    tuple val("${task.process}"), val('swave'), eval("swave-extract-alleles --version 2>/dev/null | tail -n1"), emit: versions_swave, topic: versions

    script:
    def args = task.ext.args ?: ''
    def bed_param = bed ? "--bed ${bed}" : ""
    def vcf_param = vcf ? "--vcf ${vcf}" : ""
    def is_ref_flag = is_ref ? "--is_ref" : ""
    def ref_bed_param = (ref_bed && !is_ref) ? "--ref_bed ${ref_bed}" : ""
    """
    swave-extract-alleles \\
        --gfa_fasta $gfa_fasta \\
        --sample_id ${meta.id} \\
        ${bed_param} \\
        ${vcf_param} \\
        ${is_ref_flag} \\
        ${ref_bed_param} \\
        ${args} \\
        --output_dir .
    """

    stub:
    """
    touch ${meta.id}_alleles.fa
    touch ${meta.id}_equal_paths.txt
    """
}