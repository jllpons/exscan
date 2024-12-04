process BEDTOOLS_GETFASTA {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_2@sha256:bf4b0fa21a081e5e071f273787b1e67942d331df' :
        'quay.io/biocontainers/bedtools@sha256:38756b5ac5d0368e91e85a3ed80cc40827506ebd63d449f5418befcba899b486' }"

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
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_2@sha256:bf4b0fa21a081e5e071f273787b1e67942d331df' :
        'quay.io/biocontainers/bedtools@sha256:38756b5ac5d0368e91e85a3ed80cc40827506ebd63d449f5418befcba899b486' }"

    publishDir "${params.outdir}/bedtools", mode: 'copy', overwrite: true, pattern: '*intersect.out'

    input:
    path gff_a
    path gff_b

    output:
    path 'bedtools_intersect.out', emit: bedtools_intersect
    path 'versions.yml',           emit: versions

    script:
    """
    bedtools intersect -wa -wb -a ${gff_a} -b ${gff_b} -filenames > bedtools_intersect.out

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        bedtools: \$(bedtools --version | sed -e 's/bedtools v//g')
    END_VERSIONS
    """
}

