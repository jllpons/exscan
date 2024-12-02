process GAWK_HEADER_METADATA {
    label "process_single"

    publishDir "${params.outdir}/gawk", mode: 'copy', overwrite: true, pattern: '*.translated.fasta'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'quay.io/biocontainers/gawk@sha256:701d6199235b36d054c24b1d0a889ca5e9740e301e4b46651f54d59576b73cd0' }"

    input:
    path fasta

    output:
    path '*.translated.fasta', emit: fasta_translated
    path 'versions.yml',        emit: versions

    script:
    """
    mock_seqkit_tranlsate_header_metadata.awk ${fasta} > ${fasta.baseName}.translated.fasta

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        gawk: \$(gawk --version | gawk 'NR==1{print \$3}')
    END_VERSIONS
    """
}
