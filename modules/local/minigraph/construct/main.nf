process MINIGRAPH_CONSTRUCT {
    tag "${meta.id}"
    label 'process_high_single_task'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/minigraph:0.20--he4a0461_2'
        : 'quay.io/biocontainers/minigraph:0.20--he4a0461_2'}"

    input:
    tuple val(meta), path(reference)
    path assemblies

    output:
    tuple val(meta), path("${reference.baseName.toString().split(/[_.]/)[0]}_pangenome.gfa"), emit: gfa
    path("${reference.baseName.toString().split(/[_.]/)[0]}_pangenome_info.txt"), emit: info
    tuple val("${task.process}"), val('minigraph'), eval('minigraph --version 2>&1'), emit: versions, topic: versions

    script:
    def args = task.ext.args ?: ''
    prefix = reference.baseName.toString().split(/[_.]/)[0]
    """
    minigraph \\
        -cxggs \\
        -t $task.cpus \\
        $args \\
        $reference \\
        $assemblies \\
        > ${prefix}_pangenome.gfa
    
    echo "Creation Date: \$(date)" >> ${prefix}_pangenome_info.txt
    echo "Reference Genome: $reference" >> ${prefix}_pangenome_info.txt
    echo "" >> ${prefix}_pangenome_info.txt
    echo "Involved Assemblies:" >> ${prefix}_pangenome_info.txt
    
    for file in $assemblies; do
        echo "\$file" >> ${prefix}_pangenome_info.txt
    done
    """

}