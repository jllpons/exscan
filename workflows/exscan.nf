/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HMMSCAN } from '../modules/local/hmmscan'
include { JQ_QUERY_ID } from '../modules/local/jq'
include { SEQKIT_GREP } from '../modules/local/seqkit'
include { SEQKIT_TRANSLATE } from '../modules/local/seqkit'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DOMTBLOP_DEFAULT } from '../subworkflows/domtblop'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow EXSCAN {

    take:
    ch_fasta      // str     : 'path/to/nucleotide.fasta'
    hmmdb_file    // str     : Name of HMM database file
    ch_hmmdb_dir  // str     : 'path/to/hmmdb_dir'
    ch_group      // channel : 'path/to/group'
    ch_versions   // channel : [ path(versions.yml) ]

    main:

    // DESC: Translate nucleotide sequences to amino acid sequences for all 6 reading frames.
    // ARGS: ch_fasta (str) - path to nucleotide fasta file
    // RETS: ch_translated (channel) - channel containing path to translated ORFs
    //       ch_versions (channel) - channel containing path to versions.yml
    SEQKIT_TRANSLATE(
        ch_fasta
    )
    // After each process, software versions are collected
    ch_versions = ch_versions.mix(SEQKIT_TRANSLATE.out.versions)

    // DESC: Query each ORF against a domain profile HMM database to identify
    //       potential homologous domains.
    // ARGS: hmmdb_file (str) - Name of HMM database file
    //       ch_hmmdb_dir (str) - path to HMM database directory
    //       ch_orfs_pp (channel) - channel containing path to a structured ORF file
    // RETS: hmmscan_domtblout (path) - path to hmmscan domtblout file
    HMMSCAN(
        hmmdb_file = hmmdb_file,
        ch_hmmdb_dir = ch_hmmdb_dir,
        ch_translated = SEQKIT_TRANSLATE.out.translated
    )
    ch_versions = ch_versions.mix(HMMSCAN.out.versions)


    // DESC: Parse hmmscan domtblout and perform different opeations. (See `./subworkflows/domtblop.nf`)
    // ARGS: ch_fasta (str) - path to nucleotide fasta file
    //       ch_translated (channel) - channel containing path to translated ORFs
    //       ch_domtblout (path) - path to hmmscan domtblout file
    //       ch_group (channel) - channel containing path to group file
    // RETS: hmmscan_domtblout (path) - path to hmmscan domtblout file
    DOMTBLOP_DEFAULT(
        ch_fasta = ch_fasta,
        ch_translated = SEQKIT_TRANSLATE.out.translated,
        ch_domtblout = HMMSCAN.out.hmmscan_domtblout,
        ch_group = ch_group,
        ch_versions = ch_versions
    )
    ch_versions = ch_versions.mix(DOMTBLOP_DEFAULT.out.versions)


    emit:
    hmmscan_domtblout = HMMSCAN.out.hmmscan_domtblout // path: 'path/to/hmmscan_domtblout'
    versions = ch_versions                            // channel: [ path(versions.yml) ]
}

