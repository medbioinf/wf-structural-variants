process SWAVE_WRITE_VCF {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    
    container "quay.io/swave:latest"
    
    input:
    tuple val(meta), path(tsv_files)
    path ref_fasta

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1'), emit: versions, topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 /app/swave/src/write_vcf.py \\
        --tsv_files ${tsv_files.join(' ')} \\
        --ref_fasta $ref_fasta \\
        --output_vcf "${prefix}.vcf" \\
        $args
    """
}