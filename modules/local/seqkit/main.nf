process SEQKIT_TRANSLATE {
    label "SEQKIT"

    publishDir "${params.outdir}/seqkit", mode: 'copy', overwrite: true, pattern: '*.translated.fasta'

    input:
    // `path` is an input qualifier, and `fasta` is the input name
    path fasta

    output:
    // `emit` is a way of naming the output channel
    path '*.translated.fasta', emit: translated
    path 'versions.yml',       emit: versions

    script:
    // seqkit flags:
    // -F: append frame information to sequence ID
    // -s: output individual amino acid sequences separated by stop symbol "*"
    // -f 6: translate all six frames
    // -j ${task.cpus}: use number of CPUs specified in the task (see nextflow.config)
    """
    seqkit translate -F -s -f 6 -j ${task.cpus} ${fasta} > ${fasta.baseName}.translated.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}" :
        seqkit : \$(seqkit version | sed -e "s/seqkit v//g")
    END_VERSIONS
    """
}

process SEQKIT_GREP {
    label "SEQKIT"

    input:
    path headers
    path fasta

    output:
    path 'seqkit_grep.fasta',  emit: seqkit_grep_results
    path 'versions.yml',       emit: versions

    script:
    """
    seqkit grep -f ${headers} ${fasta} > seqkit_grep.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}" :
        seqkit : \$(seqkit version | sed -e "s/seqkit v//g")
    END_VERSIONS
    """
}
