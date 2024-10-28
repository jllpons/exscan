/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GETORFS; GETORFS_POSTPROCESSING } from '../modules/local/getorfs'
include { HMMSCAN } from '../modules/local/hmmscan'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow EXSCAN {

    take:
    ch_fasta     // str    : 'path/to/nucleotide.fasta'
    ch_hmmdb     // str    : 'path/to/hmmdb/'
    ch_versions  // channel: [ path(versions.yml) ]

    main:


    // DESC: From nucleotide fasta file, extract ORFs from all 6 reading frames
    //       and translate to amino acid sequences. ORF == anything within two
    //       stop codons.
    // ARGS: ch_fasta (str) - path to nucleotide fasta file
    // OUTS: None
    // RETS: ch_orfs (channel) - channel containing path to a fasta file with all ORFs
    //       ch_versions (channel) - channel containing path to versions.yml
    GETORFS(
        ch_fasta
    )
    // After each process, software versions are collected
    ch_versions = ch_versions.mix(GETORFS.out.versions)

    // DESC: Post-process ORFs to add ORF location respect to original sequence
    //       in a structured format.
    // ARGS: ch_orfs (channel) - channel containing path to a fasta file with all ORFs
    // OUTS: None
    // RETS: ch_orfs_pp (channel) - channel containing path to a structured ORF file
    //       ch_versions (channel) - channel containing path to versions.yml
    GETORFS_POSTPROCESSING(
        ch_orfs = GETORFS.out.orfs
    )
    ch_versions = ch_versions.mix(GETORFS_POSTPROCESSING.out.versions)

    // DESC: Query each ORF against a domain profile HMM database to identify
    //       potential homologous domains.
    // ARGS: ch_hmmdb (str) - path to HMM database directory
    //       ch_orfs_pp (channel) - channel containing path to a post-processed ORF file
    // OUTS: None
    // RETS: hmmscan_domtblout (path) - path to hmmscan domtblout file
    HMMSCAN(
        ch_hmmdb = ch_hmmdb,
        ch_orfs_pp = GETORFS_POSTPROCESSING.out.orfs_pp
    )
    ch_versions = ch_versions.mix(HMMSCAN.out.versions)

    emit:
    hmmscan_domtblout = HMMSCAN.out.hmmscan_domtblout // path: 'path/to/hmmscan_domtblout'
    versions = ch_versions                            // channel: [ path(versions.yml) ]
}
