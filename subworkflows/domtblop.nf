include { BEDTOOLS_GETFASTA } from '../modules/local/bedtools'
include { DOMTBLOP_ADDSEQ } from '../modules/local/domtblop'
include { DOMTBLOP_GROUP } from '../modules/local/domtblop'
include { DOMTBLOP_PARSER } from '../modules/local/domtblop'
include { DOMTBLOP_TOBED } from '../modules/local/domtblop'
include { DOMTBLOP_TOGFF } from '../modules/local/domtblop'
include { JQ_QUERY_ID } from '../modules/local/jq'
include { SEQKIT_GREP } from '../modules/local/seqkit'


// Since it might be useful to have the aa sequences within the serialized results,
// we can extract the sequences from the fasta file based on query ids of each query result.
// We can go a step further and prefilter ( jq | seqkit grep ) the sequences
// based on the query ids to reduce the number sequences that the domtblop.py
// script has to process.
workflow DOMTBLOP_ADD_AMINOACIDSEQ_WORKFLOW {

    take:
    ch_qresults_serialized // channel : [ path(qresults.json) ]
    ch_translated         // channel : [ path(translated.fasta) ]
    ch_versions           // channel : [ path(versions.yml) ]


    main:

    // DESC: Query all of IDs from serialized JSON file and save them to a file
    // ARGS: ch_qresults_serialized (channel [path(qresults.json)]) - Channel containing path to serialized JSON file
    // RETS: ch_qresults_ids (channel [path(qresults_ids.txt)]) - Channel containing path to file containing all IDs
    //       ch_versions (channel [path(versions.yml)]) - Channel containing path to versions.yml
    JQ_QUERY_ID(
        ch_qresults_serialized = ch_qresults_serialized,
    )
    ch_versions = ch_versions.mix(JQ_QUERY_ID.out.versions)

    // DESC: Extract sequences from fasta file containing all transalted ORFs
    // ARGS: ch_query_ids (channel [path(qresults_ids.txt)]) - Channel containing path to file containing all IDs
    //       ch_translated (channel [path(translated.fasta)]) - Channel containing path to translated ORFs
    // RETS: ch_seqkit_grepseq (channel [path(seqkit_grepseq.fasta)]) - Channel containing path to file containing all sequences
    //       ch_versions (channel [path(versions.yml)]) - Channel containing path to versions.yml
    SEQKIT_GREP(
        ch_qresults_ids = JQ_QUERY_ID.out.query_ids,
        ch_translated = ch_translated,
    )
    ch_versions = ch_versions.mix(SEQKIT_GREP.out.versions)

    // DESC: Add sequences to the serialized JSON file
    // ARGS: ch_qresults_serialized (channel [path(qresults.json)]) - Channel containing path to serialized JSON file
    //       ch_seqkit_grepseq (channel [path(seqkit_grepseq.fasta)]) - Channel containing path to file containing all sequences
    // RETS: ch_qresults_serialized_addseq (channel [path(qresults_addseq.json)]) - Channel containing path to serialized JSON file with sequences
    //       ch_versions (channel [path(versions.yml)]) - Channel containing path to versions.yml
    DOMTBLOP_ADDSEQ(
        ch_qresults_serialized = ch_qresults_serialized,
        ch_qresults_id_grepseq = SEQKIT_GREP.out.seqkit_grepseq,
        type = "aminoacid"
    )
    ch_versions = ch_versions.mix(DOMTBLOP_ADDSEQ.out.versions)


    emit:
    qresults_serialized = DOMTBLOP_ADDSEQ.out.qresults_serialized
    versions = ch_versions

}


workflow DOMTBLOP_ADD_NUCLEOTIDESEQ_WORKFLOW {

    take:
    ch_qresults_serialized // channel : [ path(qresults.json) ]
    ch_fasta              // channel : [ path(nucleotide.fasta) ]
    ch_versions           // channel : [ path(versions.yml) ]


    main:

    // DESC: Generate a bed file from the serialized JSON query results
    // ARGS: ch_qresults_serialized (channel [path(qresults.json)]) - Channel containing path to serialized JSON file
    // RETS: ch_bed (channel [path(qresults.bed)]) - Channel containing path to bed file
    //       ch_versions (channel [path(versions.yml)]) - Channel containing path to versions.yml
    DOMTBLOP_TOBED(
        ch_qresults_serialized = ch_qresults_serialized,
    )
    ch_versions = ch_versions.mix(DOMTBLOP_TOBED.out.versions)


    // DESC: Get the sequences from the fasta file based on the bed file
    // ARGS: ch_bed (channel [path(qresults.bed)]) - Channel containing path to bed file
    //       ch_fasta (channel [path(nucleotide.fasta)]) - Channel containing path to nucleotide fasta file
    // RETS: ch_bedtools_getfasta (channel [path(bedtools_getfasta.fasta)]) - Channel containing path to fasta file
    //       ch_versions (channel [path(versions.yml)]) - Channel containing path to versions.yml
    BEDTOOLS_GETFASTA(
        ch_fasta = ch_fasta,
        ch_bed = DOMTBLOP_TOBED.out.qresults_bed,
    )
    ch_versions = ch_versions.mix(BEDTOOLS_GETFASTA.out.versions)


    // DESC: Add nucleotide sequences to the serialized JSON file
    // ARGS: ch_qresults_serialized (channel [path(qresults.json)]) - Channel containing path to serialized JSON file
    //       ch_bedtools_getfasta (channel [path(bedtools_getfasta.fasta)]) - Channel containing path to fasta file RETS: ch_qresults_serialized_addseq (channel [path(qresults_addseq.json)]) - Channel containing path to serialized JSON file with sequences
    //       ch_versions (channel [path(versions.yml)]) - Channel containing path to versions.yml
    DOMTBLOP_ADDSEQ(
        ch_qresults_serialized = ch_qresults_serialized,
        ch_qresults_id_grepseq = BEDTOOLS_GETFASTA.out.bedtools_getfasta,
        type = "nucleotide"
    )
    ch_versions = ch_versions.mix(DOMTBLOP_ADDSEQ.out.versions)


    emit:
    qresults_serialized = DOMTBLOP_ADDSEQ.out.qresults_serialized
    versions = ch_versions

}


// This is the defalut steps where domtblop.py is used on query results table.
// The following operations are performed in the following order:
//     1. Parse and seralize the results of each query to JSON
//     2. Add sequence data to the serialized JSON file
//     3. Group hits within specified distance from each other.
workflow DOMTBLOP_DEFAULT {

    take:
    ch_fasta      // channel : [ path(path/to/nucleotide.fasta) ]
    ch_translated // channel : [ path(path/to/translated.fasta) ]
    ch_domtblout  // channel : [ path(path/to/domtblout) ]
    ch_group      // int     : Hit maximum grouping distance
    ch_versions   // channel : [ path(versions.yml) ]

    main:

    // DESC: Parse and seralize the results of each query to JSON
    // ARGS: ch_domtblout (channel [path(domtblout)]) - Channel contianing path to hmmscan domtblout file
    // RETS: ch_qresults_serialized (channel [path(qresults.json)]) - Channel containing path to serialized JSON file
    //       ch_versions (channel [path(versions.yml)]) - Channel containing path to versions.yml
    DOMTBLOP_PARSER(
        ch_domtblout = ch_domtblout,
    )
    ch_versions = ch_versions.mix(DOMTBLOP_PARSER.out.versions)


    // DESC: Add sequence data to the serialized JSON file
    // ARGS: ch_qresults_serialized (channel [path(qresults.json)]) - Channel containing path to serialized JSON file
    //       ch_translated (channel [path(translated.fasta)]) - Channel containing path to translated ORFs
    // RETS: ch_qresults_serialized_addseq (channel [path(qresults_addseq.json)]) - Channel containing path to serialized JSON file with sequences
    //       ch_versions (channel [path(versions.yml)]) - Channel containing path to versions.yml
    DOMTBLOP_ADD_AMINOACIDSEQ_WORKFLOW(
        ch_qresults_serialized = DOMTBLOP_PARSER.out.qresults_serialized,
        ch_translated = ch_translated,
        ch_versions = ch_versions,
    )
    ch_versions = ch_versions.mix(DOMTBLOP_ADD_AMINOACIDSEQ_WORKFLOW.out.versions)


    // DESC: Add nucleotide sequences to the serialized JSON file
    // ARGS: ch_qresults_serialized (channel [path(qresults.json)]) - Channel containing path to serialized JSON file
    //       ch_fasta (channel [path(nucleotide.fasta)]) - Channel containing path to nucleotide fasta file
    // RETS: ch_qresults_serialized_addseq (channel [path(qresults_addseq.json)]) - Channel containing path to serialized JSON file with sequences
    //       ch_versions (channel [path(versions.yml)]) - Channel containing path to versions.yml
    DOMTBLOP_ADD_NUCLEOTIDESEQ_WORKFLOW(
        ch_qresults_serialized = DOMTBLOP_ADD_AMINOACIDSEQ_WORKFLOW.out.qresults_serialized,
        ch_fasta = ch_fasta,
        ch_versions = ch_versions,
    )
    ch_versions = ch_versions.mix(DOMTBLOP_ADD_NUCLEOTIDESEQ_WORKFLOW.out.versions)


    // DESC: Group hits within specified distance from each other.
    // ARGS: ch_qresults_serialized (channel [path(qresults_addseq.json)]) - Channel containing path to serialized JSON file with sequences
    //       params.group (int) - Maximum distance between hits to group
    // RETS: ch_qresults_grouped (channel [path(qresults_grouped.json)]) - Channel containing path to serialized JSON file with grouped hits
    //       ch_versions (channel [path(versions.yml)]) - Channel containing path to versions.yml
    DOMTBLOP_GROUP(
        ch_hits_serialized = DOMTBLOP_ADD_NUCLEOTIDESEQ_WORKFLOW.out.qresults_serialized,
        params.group
    )
    ch_versions = ch_versions.mix(DOMTBLOP_GROUP.out.versions)


    // DESC: Convert domain aligments from serialized JSON to GFF format
    // ARGS: ch_qresults_serialzied (channel [path(qresults_grouped.json)]) - Channel containing path to serialized JSON file with grouped hits
    // RETS: ch_gff (channel [path(qresults.gff)]) - Channel containing path to GFF file
    //       ch_versions (channel [path(versions.yml)]) - Channel containing path to versions.yml
    DOMTBLOP_TOGFF(
        DOMTBLOP_GROUP.out.qresults_serialized,
    )
    ch_versions = ch_versions.mix(DOMTBLOP_TOGFF.out.versions)


    emit:
    versions = ch_versions

}
