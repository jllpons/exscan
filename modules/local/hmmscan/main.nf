process HMMSCAN {
    label 'HMMER'

    publishDir "${params.outdir}/hmmscan", mode: 'copy', overwrite: true, pattern: '*.domtblout'

    input:
    val  hmmdb_file
    path hmmdb_dir
    path translated

    output:
    path '*.domtblout',       emit: hmmscan_domtblout
    path 'hmmscan.out',       emit: hmmscan_out
    path 'versions.yml',      emit: versions

    script:
    """
    hmmscan --domtblout ${translated.baseName}.domtblout ${hmmdb_dir}/${hmmdb_file} ${translated} > hmmscan.out

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        hmmscan: \$(hmmscan -h | head -n 2 | tail -n 1 | awk '{print \$3}')
    END_VERSIONS
    """
}

