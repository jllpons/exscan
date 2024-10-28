process WRITE_VERSIONS {
    input:
    val version_info from ch_versions.collect()

    output:
    path "versions.yml" to "${params.outdir}/pipeline_info"

    script:
    """
    cat ${version_info} > versions.yml
    """
}
