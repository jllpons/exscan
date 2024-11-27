process JQ_QUERY_ID {
    label "process_single"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/jq:1.6@sha256:924b7a56a626c51d633da4c96fa703e25e6af627' :
        'quay.io/biocontainers/jq@sha256:8ead3d706ddf3c848a4ec6e5740fb4bdcd5454a97e7c533b7417acede1de4c74' }"

    input:
    path queryresults

    output:
    path 'query_ids.txt',      emit: query_ids
    path 'versions.yml',       emit: versions

    script:
    """
    jq -r '.query_id' ${queryresults} > query_ids.txt

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        jq: \$(jq -V | sed -n 's/jq-//p')
    END_VERSIONS
    """
}
