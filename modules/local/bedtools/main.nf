process BEDTOOLS_GETFASTA {
    label 'BEDTOOLS'

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

