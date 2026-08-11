process SWAVE_WRITE_VCF {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'docker://jonahkps/panswave:0.1.0'
        : 'docker.io/jonahkps/panswave:0.1.0'}"
    
    input:
    tuple val(meta), path(tsv_files)
    path ref_fasta
    path equal_paths

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val("${task.process}"), val('swave'), eval("swave-write-vcf --version 2>/dev/null | tail -n1"), emit: versions_swave, topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def equal_paths_param = equal_paths ? "--equal_paths ${equal_paths}" : ""
    """
    swave-write-vcf \\
        --tsv_files ${tsv_files.join(' ')} \\
        --ref_fasta $ref_fasta \\
        --output_vcf "${prefix}.vcf" \\
        ${equal_paths_param} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.hap_level.vcf
    touch ${prefix}.hap_level.split.vcf
    touch ${prefix}.sample_level.vcf
    touch ${prefix}.sample_level.split.vcf
    """
}