process HMMSCAN {
    label 'HMMER'

    input:
    path hmmdb
    path orfs_pp

    output:
    path 'hmmscan.domtblout', emit: hmmscan_domtblout
    path 'hmmscan.out',       emit: hmmscan_out
    path 'versions.yml',      emit: versions

    script:
    """
    hmmscan --cpu ${task.cpus} --domtblout hmmscan.domtblout ${hmmdb}/hmmdb ${orfs_pp} > hmmscan.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        "hmmscan": \$(hmmscan -h | head -n 2 | tail -n 1 | awk '{print \$2}')
    """
}

