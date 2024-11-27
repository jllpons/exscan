process SEQKIT_TRANSLATE {
    label "process_multi"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0@sha256:4a8ac48e949803a6a118fe18614d1a93ffe64ba4' :
        'quay.io/biocontainers/seqkit@sha256:548b9b2686a311feab2d6811d113f71211280ae42b9945d51060a9049b21600e' }"

    publishDir "${params.outdir}/seqkit", mode: 'copy', overwrite: true, pattern: '*translated.fasta'


    input:
    path fasta

    output:
    path '*translated.fasta',     emit: fasta_translated
    path 'versions.yml',          emit: versions

    script:
    """
    seqkit translate --append-frame --out-subseqs --frame 6 --min-len 12 --threads ${task.cpus} ${fasta} > ${fasta.baseName}.translated.fasta

    seqkit split2 --quiet --by-size 300000 --threads ${task.cpus} --by-size-prefix translated.split.fasta ${fasta.baseName}.translated.fasta --out-dir ${params.outdir}/seqkit/split

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        seqkit: \$(seqkit version | sed -e "s/seqkit v//g")
    END_VERSIONS
    """
}


process SEQKIT_SPLIT {
    label "process_multi"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0@sha256:4a8ac48e949803a6a118fe18614d1a93ffe64ba4' :
        'quay.io/biocontainers/seqkit@sha256:548b9b2686a311feab2d6811d113f71211280ae42b9945d51060a9049b21600e' }"

    publishDir "${params.outdir}/seqkit", mode: 'copy', overwrite: true, pattern: 'seqkit_split/*fasta'


    input:
    path fasta

    output:
    path 'seqkit_split/*fasta', emit: fasta_split
    path 'versions.yml',                                          emit: versions

    script:
    """
    seqkit split2 --quiet --by-size 300000 --threads ${task.cpus} --by-size-prefix translated.split.fasta ${fasta} --out-dir seqkit_split

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        seqkit: \$(seqkit version | sed -e "s/seqkit v//g")
    END_VERSIONS
    """
}


process SEQKIT_GREP {
    label "process_single"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0@sha256:4a8ac48e949803a6a118fe18614d1a93ffe64ba4' :
        'quay.io/biocontainers/seqkit@sha256:548b9b2686a311feab2d6811d113f71211280ae42b9945d51060a9049b21600e' }"


    input:
    path headers
    path fasta

    output:
    path 'seqkit_grepseq.fasta',  emit: seqkit_grepseq
    path 'versions.yml',       emit: versions

    script:
    """
    seqkit grep -f ${headers} ${fasta} > seqkit_grepseq.fasta

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        seqkit: \$(seqkit version | sed -e "s/seqkit v//g")
    END_VERSIONS
    """
}
