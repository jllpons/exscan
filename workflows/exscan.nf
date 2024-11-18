/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DOMTBLOP_ADDSEQ } from '../modules/local/domtblop'
include { DOMTBLOP_PARSER } from '../modules/local/domtblop'
include { DOMTBLOP_GROUP } from '../modules/local/domtblop'
include { HMMSCAN } from '../modules/local/hmmscan'
include { JQ_QUERY_ID } from '../modules/local/jq'
include { SEQKIT_GREP } from '../modules/local/seqkit'
include { SEQKIT_TRANSLATE } from '../modules/local/seqkit'

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
    // OUTS: None
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
    // OUTS: None
    // RETS: hmmscan_domtblout (path) - path to hmmscan domtblout file
    HMMSCAN(
        hmmdb_file = hmmdb_file,
        ch_hmmdb_dir = ch_hmmdb_dir,
        ch_translated = SEQKIT_TRANSLATE.out.translated
    )
    ch_versions = ch_versions.mix(HMMSCAN.out.versions)

    // DESC: Parse and serialize to JSON the hmmscan domtblout file.
    // ARGS: ch_domtblout (channel) - channel containing path to hmmscan domtblout file
    // OUTS: None
    // RETS: ch_hits_serialized (channel) - channel containing path to serialized hits
    DOMTBLOP_PARSER(
        ch_domtblout = HMMSCAN.out.hmmscan_domtblout
    )
    ch_versions = ch_versions.mix(DOMTBLOP_PARSER.out.versions)


    // NOTE: Next two steps are not strictly necessary.
    //       Since it might be useful to have the sequences within the serialized results,
    //       we can extract the sequences from the fasta file based on the query ids.
    //       We go a step further and prefilter the sequences based on the query ids
    //       to reduce the number of sequences to be handled.


    // DESC: Query the query ids of hmmscan results and save it to a list.
    // ARGS: ch_hits_serialized (channel) - channel containing path to serialized hits
    // OUTS: None
    // RETS: ch_query_ids (channel)
    JQ_QUERY_ID(
        ch_hits_serialized = DOMTBLOP_PARSER.out.queryresults
    )
    ch_versions = ch_versions.mix(JQ_QUERY_ID.out.versions)

    // DESC: Extract sequences from the fasta file based on the query ids.
    // ARGS: ch_query_ids (channel) - channel containing query ids
    //       ch_translated (channel) - channel containing path to translated ORFs
    // OUTS: None
    // RETS: ch_hits_fasta (channel) - channel containing path to extracted sequences
    SEQKIT_GREP(
        ch_query_ids = JQ_QUERY_ID.out.query_ids,
        ch_fasta = SEQKIT_TRANSLATE.out.translated
    )
    ch_versions = ch_versions.mix(SEQKIT_GREP.out.versions)

    // DESC: Add sequences to the serialized hits.
    // ARGS: ch_hits_serialized (channel) - channel containing path to serialized hits
    //       ch_hits_fasta (channel) - channel containing path to extracted sequences
    // OUTS: None
    // RETS: ch_hits_serialized (channel) - channel containing path to serialized hits
    DOMTBLOP_ADDSEQ(
        ch_hits_serialized = DOMTBLOP_PARSER.out.queryresults,
        ch_hits_fasta = SEQKIT_GREP.out.seqkit_grep_results
    )
    ch_versions = ch_versions.mix(DOMTBLOP_ADDSEQ.out.versions)


    if(params.group != null) {
        // DESC: Group hits within specified distance from each other.
        // ARGS: ch_hits_serialized (channel) - channel containing path to serialized hits
        // OUTS: None
        // RETS: ch_hits_grouped (channel) - channel containing path to grouped hits
        DOMTBLOP_GROUP(
            ch_hits_serialized = DOMTBLOP_PARSER.out.queryresults,
            params.group
        )
        ch_versions = ch_versions.mix(DOMTBLOP_GROUP.out.versions)
    }

    emit:
    hmmscan_domtblout = HMMSCAN.out.hmmscan_domtblout // path: 'path/to/hmmscan_domtblout'
    versions = ch_versions                            // channel: [ path(versions.yml) ]
}

