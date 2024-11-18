process DOMTBLOP_PARSER {
    label 'BIOPYTHON'

    publishDir "${params.outdir}/domtblop", mode: 'copy', overwrite: true, pattern: "*queryresults.json"

    input:
    path domtblout

    output:
    path '*queryresults.json',  emit: queryresults
    path 'versions.yml',        emit: versions

    script:
    """
    domtblop.py parser ${domtblout} > ${domtblout.baseName}.queryresults.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        "Python": \$(python -V | awk '{print \$2}')
        "Biopython": \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}


process DOMTBLOP_GROUP {
    label 'BIOPYTHON'

    publishDir "${params.outdir}/domtblop", mode: 'copy', overwrite: true, pattern: "*grouped.json"

    input:
    path queryresults
    val group

    output:
    path '*grouped.json',       emit: queryresults
    path 'versions.yml',        emit: versions

    script:
    """
    domtblop.py group -d ${group} ${queryresults} > ${queryresults.baseName}.grouped.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        "Python": \$(python -V | awk '{print \$2}')
        "Biopython": \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}

process DOMTBLOP_ADDSEQ {
    label 'BIOPYTHON'

    publishDir "${params.outdir}/domtblop", mode: 'copy', overwrite: true, pattern: "*addseq.json"

    input:
    path queryresults
    path fasta

    output:
    path '*addseq.json',        emit: queryresults
    path 'versions.yml',        emit: versions

    script:
    """
    domtblop.py addseq ${queryresults} ${fasta} > ${queryresults.baseName}.addseq.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        "Python": \$(python -V | awk '{print \$2}')
        "Biopython": \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
