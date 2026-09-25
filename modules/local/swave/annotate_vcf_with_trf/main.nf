process SWAVE_ANNOTATE_VCF_WITH_TRF {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'docker://jonahkps/panswave:0.1.0'
        : 'docker.io/jonahkps/panswave:0.1.0'}"

    input:
    tuple val(meta), path(vcf)
    path(dat_files)
    path(fasta_files)

    output:
    tuple val(meta), path("*.trf_annotated.vcf"), emit: vcf
    tuple val("${task.process}"), val('swave'), eval("swave-annotate-vcf-with-trf --version 2>/dev/null | tail -n1"), emit: versions_swave, topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    swave-annotate_vcf_with_trf \\
        --vcf $vcf \\
        --dat $dat_files \\
        --fasta $fasta_files \\
        --output_vcf ${prefix}.trf_annotated.vcf \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.trf_annotated.vcf
    """
}