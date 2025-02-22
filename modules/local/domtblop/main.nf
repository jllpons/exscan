process DOMTBLOP_ADDSEQ {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.79' :
        'quay.io/biocontainers/biopython:1.79' }"

    publishDir "${params.outdir}/domtblop", mode: 'copy', overwrite: true, pattern: "*addseq.json"

    input:
    path qresults_serialized
    path fasta
    val  type

    output:
    path '*addseq.json',        emit: qresults_serialized
    path 'versions.yml',        emit: versions

    script:
    if (type == 'nucleotide')
        """
        domtblop.py addseq ${qresults_serialized} ${fasta} --seq-type ${type} > ${qresults_serialized.baseName}.addseq.json

        cat <<-END_VERSIONS > versions.yml
        ${task.process}:
            Python: \$(python -V | awk '{print \$2}')
            Biopython: \$(python -c "import Bio; print(Bio.__version__)")
        END_VERSIONS
        """
    else
        """
        domtblop.py addseq ${qresults_serialized} ${fasta} --seq-type ${type} > ${qresults_serialized.baseName}.addseq.json

        cat <<-END_VERSIONS > versions.yml
        ${task.process}:
            Python: \$(python -V | awk '{print \$2}')
            Biopython: \$(python -c "import Bio; print(Bio.__version__)")
        END_VERSIONS
        """
}

// WARNING: Processes in the workflow cannot share the same name, so if used,
//          different filtering processes should be encapsulated in different subworkflows.
//          Not used atm
process DOMTBLOP_FILTER {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.79' :
        'quay.io/biocontainers/biopython:1.79' }"

    publishDir "${params.outdir}/domtblop", mode: 'copy', overwrite: true, pattern: "*.filter.json"

    input:
    path qresults_serialized
    val  filter
    val  value

    output:
    path '*.filter.json', emit: qresults_serialized
    path 'versions.yml', emit: versions

    script:

    switch(filter){
        case 'dom-ievalue':
            """
            domtblop.py filter --dom-ievalue ${value} ${qresults_serialized} > ${qresults_serialized.baseName}.filter.json
            """

        case 'min-alignment-length':
            """
            domtblop.py filter --min-alignment-length ${value} ${qresults_serialized} > ${qresults_serialized.baseName}.filter.json
            """

        case 'keep-only-intersect':
            """
            domtblop.py filter --keep-only-intersect ${qresults_serialized} > ${qresults_serialized.baseName}.filter.json
            """

        default:
            log.error("Invalid filter: ${filter}")
    }
}


process DOMTBLOP_FILTER_BY_DOMAIN_IEVALUE {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.79' :
        'quay.io/biocontainers/biopython:1.79' }"

    publishDir "${params.outdir}/domtblop", mode: 'copy', overwrite: true, pattern: "*.filterDomIEvalue.json"

    input:
    path qresults_serialized
    val threshold
    val keep_best_hit

    output:
    path '*.filterDomIEvalue.json', emit: qresults_serialized
    path 'versions.yml',        emit: versions

    script:

    if (keep_best_hit)
        """
        domtblop.py filter --best-hit --dom-ievalue ${threshold} ${qresults_serialized} \
            | domtblop.py filter --best-hit --seq-evalue "1.0" > ${qresults_serialized.baseName}.filterDomIEvalue.json
        cat <<-END_VERSIONS > versions.yml
        ${task.process}:
            Python: \$(python -V | awk '{print \$2}')
            Biopython: \$(python -c "import Bio; print(Bio.__version__)")
        END_VERSIONS
        """
    else
        """
        domtblop.py filter --dom-ievalue ${threshold} ${qresults_serialized} > ${qresults_serialized.baseName}.filterDomIEvalue.json

        cat <<-END_VERSIONS > versions.yml
        ${task.process}:
            Python: \$(python -V | awk '{print \$2}')
            Biopython: \$(python -c "import Bio; print(Bio.__version__)")
        END_VERSIONS
        """
}


process DOMTBLOP_FILTER_BY_MIN_ALIGNMENT_LENGTH {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.79' :
        'quay.io/biocontainers/biopython:1.79' }"

    publishDir "${params.outdir}/domtblop", mode: 'copy', overwrite: true, pattern: "*filterMinAlignmentLength.json"

    input:
    path qresults_serialized
    val  min_alignment_length

    output:
    path '*filterMinAlignmentLength.json', emit: qresults_serialized
    path 'versions.yml',                   emit: versions

    script:
    """
    domtblop.py filter --min-alignment-length ${min_alignment_length} ${qresults_serialized} > ${qresults_serialized.baseName}.filterMinAlignmentLength.json

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        Python: \$(python -V | awk '{print \$2}')
        Biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}

process DOMTBLOP_FILTER_KEEP_ONLY_INTERSECT {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.79' :
        'quay.io/biocontainers/biopython:1.79' }"

    publishDir "${params.outdir}/domtblop", mode: 'copy', overwrite: true, pattern: "*filterKeepOnlyIntersect.json"

    input:
    path qresults_serialized

    output:
    path '*filterMinAlignmentLength.json', emit: qresults_serialized
    path 'versions.yml',                   emit: versions

    script:
    """
    domtblop.py filter --keep-only-intersect ${qresults_serialized} > ${qresults_serialized.baseName}.filterKeepOnlyIntersect.json

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        Python: \$(python -V | awk '{print \$2}')
        Biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}


process DOMTBLOP_GFFINTERSECT {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.79' :
        'quay.io/biocontainers/biopython:1.79' }"

    publishDir "${params.outdir}/domtblop", mode: 'copy', overwrite: true, pattern: "*gffintersect.json"

    input:
    path qresults_serialized
    path bedtools_intersect_out

    output:
    path '*gffintersect.json',  emit: qresults_serialized
    path 'versions.yml',        emit: versions

    script:
    """
    domtblop.py gffintersect ${qresults_serialized} ${bedtools_intersect_out} > ${qresults_serialized.baseName}.gffintersect.json

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        Python: \$(python -V | awk '{print \$2}')
        Biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}


process DOMTBLOP_GROUP {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.79' :
        'quay.io/biocontainers/biopython:1.79' }"

    publishDir "${params.outdir}/domtblop", mode: 'copy', overwrite: true, pattern: "*group.json"

    input:
    path in_qresults_serialized
    val group

    output:
    path '*group.json',         emit: qresults_serialized
    path 'versions.yml',        emit: versions

    script:
    """
    domtblop.py group -d ${group} ${in_qresults_serialized} > ${in_qresults_serialized.baseName}.group.json

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        Python: \$(python -V | awk '{print \$2}')
        Biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}


process DOMTBLOP_PARSER {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.79' :
        'quay.io/biocontainers/biopython:1.79' }"

    publishDir "${params.outdir}/domtblop", mode: 'copy', overwrite: true, pattern: "*qresults_serialized.json"

    input:
    path hmmscan_domtblout
    val  sequence_type

    output:
    path '*qresults_serialized.json',  emit: qresults_serialized
    path 'versions.yml',        emit: versions

    script:
    """
    domtblop.py parser ${hmmscan_domtblout} --seq-type ${sequence_type} > ${hmmscan_domtblout.baseName}.qresults_serialized.json

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        Python: \$(python -V | awk '{print \$2}')
        Biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}


process DOMTBLOP_TOBED {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.79' :
        'quay.io/biocontainers/biopython:1.79' }"

    publishDir "${params.outdir}/domtblop", mode: 'copy', overwrite: true, pattern: "*bed"

    input:
    path qresults_serialized

    output:
    path '*bed',                emit: qresults_bed
    path 'versions.yml',        emit: versions

    script:
    """
    domtblop.py tobed ${qresults_serialized} > ${qresults_serialized.baseName}.bed

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        Python: \$(python -V | awk '{print \$2}')
        Biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}


process DOMTBLOP_TOCSV {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.79' :
        'quay.io/biocontainers/biopython:1.79' }"

    publishDir "${params.outdir}/domtblop", mode: 'copy', overwrite: true, pattern: "*.csv"

    input:
    path qresults_serialized

    output:
    path '*sequenceCentric.csv', emit: qresults_sequence_csv
    path '*domainCentric.csv',   emit: qresults_domain_csv
    //path '*groupCentric.csv',    emit: qresults_group_csv
    //path '*gffIntersect.csv',    emit: qresults_gff_csv
    //path '*intersectionMap.csv', emit: qresults_map_csv
    path 'versions.yml',        emit: versions

    script:
    """
    domtblop.py tocsv ${qresults_serialized} --kind sequence > ${qresults_serialized.baseName}.sequenceCentric.csv
    domtblop.py tocsv ${qresults_serialized} --kind domain   > ${qresults_serialized.baseName}.domainCentric.csv


    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        Python: \$(python -V | awk '{print \$2}')
        Biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
    /*
    domtblop.py tocsv ${qresults_serialized} --kind group    > ${qresults_serialized.baseName}.groupCentric.csv
    domtblop.py tocsv ${qresults_serialized} --kind gff      > ${qresults_serialized.baseName}.gffIntersect.csv
    domtblop.py tocsv ${qresults_serialized} --kind map      > ${qresults_serialized.baseName}.intersectionMap.csv
    */
}


process DOMTBLOP_TOFASTA {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.79' :
        'quay.io/biocontainers/biopython:1.79' }"

    publishDir "${params.outdir}/domtblop", mode: 'copy', overwrite: true, pattern: "*.fasta"

    input:
    path qresults_serialized

    output:
    path '*.fasta',             emit: qresults_fasta
    path 'versions.yml',        emit: versions

    script:
    """
    domtblop.py tofasta ${qresults_serialized} --seq-kind nucleotide > ${qresults_serialized.baseName}.nucleotide.fasta
    domtblop.py tofasta ${qresults_serialized} --seq-kind protein > ${qresults_serialized.baseName}.protein.fasta
    domtblop.py tofasta ${qresults_serialized} --seq-kind domain-alignment > ${qresults_serialized.baseName}.domain_alignment.fasta

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        Python: \$(python -V | awk '{print \$2}')
        Biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}


process DOMTBLOP_TOGFF {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.79' :
        'quay.io/biocontainers/biopython:1.79' }"

    publishDir "${params.outdir}/domtblop", mode: 'copy', overwrite: true, pattern: "*gff"

    input:
    path qresults_serialized

    output:
    path '*gff',                emit: qresults_gff
    path 'versions.yml',        emit: versions

    script:
    """
    domtblop.py togff ${qresults_serialized} > ${qresults_serialized.baseName}.gff

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        Python: \$(python -V | awk '{print \$2}')
        Biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}

