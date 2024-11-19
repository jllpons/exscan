process JQ_QUERY_ID {
    label "JQ"

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
