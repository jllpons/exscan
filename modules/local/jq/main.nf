process JQ_QUERY_ID {
    label "process_single"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/jq:1.6' :
        'quay.io/biocontainers/jq:1.6' }"

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
