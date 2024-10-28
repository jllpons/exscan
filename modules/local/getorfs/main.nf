process GETORFS {
    label 'EMBOSS'

    input:
    // `path` is an input qualifier, and `fasta` is the input name
    path fasta

    output:
    // `emit` is a way of naming the output channel
    path '*.orfs.fasta', emit: orfs
    path 'versions.yml', emit: versions

    script:
    """
    getorf \\
        -auto \\
        -sequence ${fasta} \\
        -outseq ${fasta.baseName}.orfs.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getorf: \$(getorf -version 2>&1 | sed '/EMBOSS:/!d; s/.*://')
    END_VERSIONS
    """
}

process GETORFS_POSTPROCESSING {
    label 'GAWK'

    input:
    // `orfs` is the output channel of the previous process
    path orfs

    output:
    path '*.pp.fasta', emit: orfs_pp
    path 'versions.yml', emit: versions

    script:
    """
    getorf_pp.awk $orfs > ${orfs.baseName}.pp.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(gawk --version | head -n 1 | gawk '{print \$3}' | sed 's/,//')
    END_VERSIONS
    """
}
