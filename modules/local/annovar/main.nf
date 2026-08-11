process ANNOVAR {
    tag "${meta.id}"
    label 'process_medium'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'docker://docker.io/perl:5.32'
        : 'docker.io/perl:5.32'}"

    input:
    tuple val(meta), path(vcf)
    path annovar_dir
    val buildver

    output:
    tuple val(meta), path("*multianno.vcf"), emit: vcf
    tuple val("${task.process}"), val('annovar'), eval("stat -c %y ${annovar_dir}/table_annovar.pl | cut -d' ' -f1"), emit: versions_annovar, topic: versions

script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    perl ${annovar_dir}/table_annovar.pl \\
        ${vcf} \\
        ${annovar_dir}/humandb/ \\
        -buildver ${buildver} \\
        -out ${prefix} \\
        -vcfinput \\
        -protocol refGene \\
        -operation g \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.${buildver}_multianno.vcf
    """
}