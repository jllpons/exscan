process DOMTBLOP_PARSER {
    label 'BIOPYTHON'

    publishDir "${params.outdir}/domtblop", mode: 'copy', overwrite: true, pattern: "*queryresults.json"

    input:
    path domtblout

    output:
    path '*queryresults.json',  emit: qresults_serialized
    path 'versions.yml',        emit: versions

    script:
    """
    domtblop.py parser ${domtblout} > ${domtblout.baseName}.queryresults.json

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        Python: \$(python -V | awk '{print \$2}')
        Biopython: \$(python -c "import Bio; print(Bio.__version__)")
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
    path '*grouped.json',       emit: qresults_serialized
    path 'versions.yml',        emit: versions

    script:
    """
    domtblop.py group -d ${group} ${queryresults} > ${queryresults.baseName}.grouped.json

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        Python: \$(python -V | awk '{print \$2}')
        Biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}

process DOMTBLOP_ADDSEQ {
    label 'BIOPYTHON'

    publishDir "${params.outdir}/domtblop", mode: 'copy', overwrite: true, pattern: "*addseq.json"

    input:
    path queryresults
    path fasta
    val  type

    output:
    path '*addseq.json',        emit: qresults_serialized
    path 'versions.yml',        emit: versions

    script:
    if (type == 'nucleotide')
        """
        domtblop.py addseq --nucleotide ${queryresults} ${fasta} > ${queryresults.baseName}.addseq.json

        cat <<-END_VERSIONS > versions.yml
        ${task.process}:
            Python: \$(python -V | awk '{print \$2}')
            Biopython: \$(python -c "import Bio; print(Bio.__version__)")
        END_VERSIONS
        """
    else
        """
        domtblop.py addseq ${queryresults} ${fasta} > ${queryresults.baseName}.addseq.json

        cat <<-END_VERSIONS > versions.yml
        ${task.process}:
            Python: \$(python -V | awk '{print \$2}')
            Biopython: \$(python -c "import Bio; print(Bio.__version__)")
        END_VERSIONS
        """
}

process DOMTBLOP_TOGFF {
    label 'BIOPYTHON'

    publishDir "${params.outdir}/domtblop", mode: 'copy', overwrite: true, pattern: "*gff"

    input:
    path queryresults

    output:
    path '*gff',                emit: qresults_gff
    path 'versions.yml',        emit: versions

    script:
    """
    domtblop.py togff ${queryresults} > ${queryresults.baseName}.gff

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        Python: \$(python -V | awk '{print \$2}')
        Biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}

process DOMTBLOP_TOBED {
    label 'BIOPYTHON'

    publishDir "${params.outdir}/domtblop", mode: 'copy', overwrite: true, pattern: "*bed"

    input:
    path queryresults

    output:
    path '*bed',                emit: qresults_bed
    path 'versions.yml',        emit: versions

    script:
    """
    domtblop.py tobed ${queryresults} > ${queryresults.baseName}.bed

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        Python: \$(python -V | awk '{print \$2}')
        Biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
