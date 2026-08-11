process GFA_SET_REFERENCE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
                ? 'https://depot.galaxyproject.org/singularity/gawk:5.3.0'
                : 'biocontainers/gawk:5.3.0'}"

    input:
    tuple val(meta), path(gfa)
    path(ref_contig_names)

    output:
    tuple val(meta), path("*.reftagged.gfa"), emit: gfa
    tuple val("${task.process}"), val('gawk'), eval("awk --version 2>&1 | head -n1 | sed 's/GNU Awk //; s/,.*//'"), emit: versions_gawk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    RS_TAGS=\$(cut -d'#' -f1 ${ref_contig_names} | sort -u | awk '{printf "RS:Z:%s\\t", \$1}' | sed 's/\\t\$//')

    awk -v OFS='\\t' -v rs="\$RS_TAGS" '
    NR==1 && \$1=="H" {
        if (\$0 ~ /RS:Z:/) { print; next }
        print \$0, rs; next
    }
    FNR==1 && \$1!="H" { print "H", "VN:Z:1.0", rs; print; next }
    { print }
    ' ${gfa} > ${prefix}.reftagged.gfa
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.reftagged.gfa
    """
}