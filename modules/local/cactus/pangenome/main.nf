process CACTUS_PANGENOME {
    tag "${meta.id}"
    label 'process_high'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'docker://quay.io/comparative-genomics-toolkit/cactus:v2.9.3'
        : 'quay.io/comparative-genomics-toolkit/cactus:v2.9.3'}"

    input:
    tuple val(meta), path(seqfile)
    path(fastas)

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.gfa.gz"), emit: gfa
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.sv.gfa.gz"), emit: sv_gfa, optional: true
    tuple val("${task.process}"), val('cactus'), eval("cactus --version 2>&1 | head -n1"), emit: versions_cactus, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export HOME=\$PWD
    mkdir -p ./cactus_tmp

    cactus-pangenome \\
        ./jobstore \\
        ${seqfile} \\
        --outDir . \\
        --outName ${prefix} \\
        --reference ${meta.sample} \\
        --gfa \\
        --workDir ./cactus_tmp \\
        --maxCores ${task.cpus} \\
        --maxMemory ${task.memory.toGiga()}G \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export HOME=\$PWD
    echo "" | gzip -c > ${prefix}.gfa.gz
    """
}