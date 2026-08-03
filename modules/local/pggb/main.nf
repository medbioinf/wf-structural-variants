process PGGB {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pggb:0.7.4--h9ee0642_0'
        : 'quay.io/biocontainers/pggb:0.7.4--h9ee0642_0'}"
    
    input:
    tuple val(meta), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.gfa.gz"), emit: gfa
    tuple val(meta), path("*.og"), emit: og, optional: true
    tuple val(meta), path("*.vcf.gz"), emit: vcf, optional: true
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi, optional: true
    tuple val("${task.process}"), val('pggb'), eval("pggb --version 2>&1 | head -n 1 | sed 's/pggb //g'"), emit: versions_pggb, topic: versions

    script:
    def args = task.ext.args ?: ''
    """
    pggb \\
        -i ${fasta} \\
        -o . \\
        -t ${task.cpus} \\
        ${args}

    bgzip -f *.gfa
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | bgzip -c > ${prefix}.gfa.gz
    touch ${prefix}.og
    echo "" | bgzip -c > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}