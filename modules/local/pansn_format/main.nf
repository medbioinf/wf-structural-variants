process PANSN_FORMAT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/gawk:5.3.0'
        : 'biocontainers/gawk:5.3.0'}"
    
    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.pansn.fa"), emit: fasta
    tuple val("${task.process}"), val('gawk'), eval("awk --version 2>&1 | head -n1 | sed 's/GNU Awk //; s/,.*//'"), emit: versions_gawk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def sample = meta.sample ?: meta.id
    def haplotype = meta.haplotype ?: "0"
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [[ "${fasta}" == *.gz ]]; then
        READ_CMD="gzip -cd"
    else
        READ_CMD="cat"
    fi

    \$READ_CMD ${fasta} | awk -v sample="${sample}" -v hap="${haplotype}" '
        /^>/ {
            header = \$1
            sub(/^>/, "", header)
            n = split(header, a, "#")
            contig = a[n]
            print ">" sample "#" hap "#" contig
            next
        }
        { print }
    ' > ${prefix}.pansn.fa
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pansn.fa
    """
}