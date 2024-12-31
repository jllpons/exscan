process DOMTBLOP_ADDSEQ {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.79@sha256:dc432d0b398037b797d6981ec338522e5417bbf4' :
        'quay.io/biocontainers/biopython@sha256:937556be7fd782859ece3138e0b8beae3f4645ae8c8fcf304bd56d06084ae37b' }"

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
        'https://depot.galaxyproject.org/singularity/biopython:1.79@sha256:dc432d0b398037b797d6981ec338522e5417bbf4' :
        'quay.io/biocontainers/biopython@sha256:937556be7fd782859ece3138e0b8beae3f4645ae8c8fcf304bd56d06084ae37b' }"

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
        'https://depot.galaxyproject.org/singularity/biopython:1.79@sha256:dc432d0b398037b797d6981ec338522e5417bbf4' :
        'quay.io/biocontainers/biopython@sha256:937556be7fd782859ece3138e0b8beae3f4645ae8c8fcf304bd56d06084ae37b' }"

    publishDir "${params.outdir}/domtblop", mode: 'copy', overwrite: true, pattern: "*.filterDomIEvalue.json"

    input:
    path qresults_serialized
    val threshold

    output:
    path '*.filterDomIEvalue.json', emit: qresults_serialized
    path 'versions.yml',        emit: versions

    script:
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
        'https://depot.galaxyproject.org/singularity/biopython:1.79@sha256:dc432d0b398037b797d6981ec338522e5417bbf4' :
        'quay.io/biocontainers/biopython@sha256:937556be7fd782859ece3138e0b8beae3f4645ae8c8fcf304bd56d06084ae37b' }"

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
        'https://depot.galaxyproject.org/singularity/biopython:1.79@sha256:dc432d0b398037b797d6981ec338522e5417bbf4' :
        'quay.io/biocontainers/biopython@sha256:937556be7fd782859ece3138e0b8beae3f4645ae8c8fcf304bd56d06084ae37b' }"

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
        'https://depot.galaxyproject.org/singularity/biopython:1.79@sha256:dc432d0b398037b797d6981ec338522e5417bbf4' :
        'quay.io/biocontainers/biopython@sha256:937556be7fd782859ece3138e0b8beae3f4645ae8c8fcf304bd56d06084ae37b' }"

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
        'https://depot.galaxyproject.org/singularity/biopython:1.79@sha256:dc432d0b398037b797d6981ec338522e5417bbf4' :
        'quay.io/biocontainers/biopython@sha256:937556be7fd782859ece3138e0b8beae3f4645ae8c8fcf304bd56d06084ae37b' }"

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
        'https://depot.galaxyproject.org/singularity/biopython:1.79@sha256:dc432d0b398037b797d6981ec338522e5417bbf4' :
        'quay.io/biocontainers/biopython@sha256:937556be7fd782859ece3138e0b8beae3f4645ae8c8fcf304bd56d06084ae37b' }"

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
        'https://depot.galaxyproject.org/singularity/biopython:1.79@sha256:dc432d0b398037b797d6981ec338522e5417bbf4' :
        'quay.io/biocontainers/biopython@sha256:937556be7fd782859ece3138e0b8beae3f4645ae8c8fcf304bd56d06084ae37b' }"

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


process DOMTBLOP_TOFASTA {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.79@sha256:dc432d0b398037b797d6981ec338522e5417bbf4' :
        'quay.io/biocontainers/biopython@sha256:937556be7fd782859ece3138e0b8beae3f4645ae8c8fcf304bd56d06084ae37b' }"

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
        'https://depot.galaxyproject.org/singularity/biopython:1.79@sha256:dc432d0b398037b797d6981ec338522e5417bbf4' :
        'quay.io/biocontainers/biopython@sha256:937556be7fd782859ece3138e0b8beae3f4645ae8c8fcf304bd56d06084ae37b' }"

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

