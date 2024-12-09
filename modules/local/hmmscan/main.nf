process HMMSCAN {
    label 'process_single_more_memory'


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.4--hdbdd923_2@sha256:145286a01c4352ca95d3ff76072ec37fafd9a86e' :
        'quay.io/biocontainers/hmmer@sha256:f482523ecef8ddb9bf2d63e896b03173a51e9698d6f4e5c8f03c8aefedcaedc6' }"

    publishDir "${params.outdir}/hmmscan", mode: 'copy', overwrite: true, pattern: '*.domtblout'


    input:
    val  hmmdb_file
    path hmmdb_dir
    each path(fasta)

    output:
    path '*.domtblout',       emit: hmmscan_domtblout
    path 'hmmscan.out',       emit: hmmscan_out
    path 'versions.yml',      emit: versions

    script:
    """
    hmmscan --domtblout ${fasta.baseName}.domtblout ${hmmdb_dir}/${hmmdb_file} ${fasta} > hmmscan.out

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        hmmscan: \$(hmmscan -h | head -n 2 | tail -n 1 | awk '{print \$3}')
    END_VERSIONS
    """
}


process MERGE_DOMTBLOUT_RESULTS {
    label "process_single"

    publishDir "${params.outdir}/hmmscan", mode: 'copy', overwrite: true, pattern: 'hmmmscan_merged.domtblout'


    input:
    path domtblout_files

    output:
    path 'hmmscan_merged.domtblout', emit: hmmscan_merged_domtblout

    script:
    """
    cat ${domtblout_files} | grep -v '^#' > hmmscan_merged.domtblout
    """
}

