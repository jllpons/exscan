process BEDTOOLS_GETFASTA {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--h13024bc_3' :
        'quay.io/biocontainers/bedtools:2.31.1--h13024bc_3' }"

    publishDir "${params.outdir}/bedtools", mode: 'copy', overwrite: true, pattern: '*.fasta'

    input:
    path fasta
    path bed

    output:
    path 'bedtools_getfasta.fasta',  emit: bedtools_getfasta
    path 'versions.yml',             emit: versions

    script:
    """
    bedtools getfasta -fi ${fasta} -bed ${bed} -fo bedtools_getfasta.fasta -nameOnly

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        bedtools: \$(bedtools --version | sed -e 's/bedtools v//g')
    END_VERSIONS
    """
}


process BEDTOOLS_INTERSECT {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--h13024bc_3' :
        'quay.io/biocontainers/bedtools:2.31.1--h13024bc_3' }"

    publishDir "${params.outdir}/bedtools", mode: 'copy', overwrite: true, pattern: '*intersect.out'

    input:
    path gff_a
    path gff_b

    output:
    path 'bedtools_intersect.out', emit: bedtools_intersect
    path 'versions.yml',           emit: versions

    script:
    """
    bedtools intersect -wa -wb -a ${gff_a} -b ${gff_b} > bedtools_intersect.out

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        bedtools: \$(bedtools --version | sed -e 's/bedtools v//g')
    END_VERSIONS
    """
}

