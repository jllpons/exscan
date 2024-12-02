include { BEDTOOLS_GETFASTA                       } from '../modules/local/bedtools'
include { DOMTBLOP_ADDSEQ                         } from '../modules/local/domtblop'
include { DOMTBLOP_FILTER_BY_DOMAIN_IEVALUE       } from '../modules/local/domtblop'
include { DOMTBLOP_FILTER_BY_MIN_ALIGNMENT_LENGTH } from '../modules/local/domtblop'
include { DOMTBLOP_GROUP                          } from '../modules/local/domtblop'
include { DOMTBLOP_PARSER                         } from '../modules/local/domtblop'
include { DOMTBLOP_TOBED                          } from '../modules/local/domtblop'
include { DOMTBLOP_TOFASTA                        } from '../modules/local/domtblop'
include { DOMTBLOP_TOGFF                          } from '../modules/local/domtblop'
include { JQ_QUERY_ID                             } from '../modules/local/jq'
include { SEQKIT_GREP                             } from '../modules/local/seqkit'


// Since it might be useful to have the aa sequences within the serialized results,
// we can extract the sequences from the fasta file based on query ids of each query result.
// We can go a step further and prefilter ( jq | seqkit grep ) the sequences
// based on the query ids to reduce the number sequences that the domtblop.py
// script has to process.
workflow DOMTBLOP_ADD_AMINOACIDSEQ_WORKFLOW {

    take:
    ch_qresults_serialized // channel : [ path(qresults.json) ]
    ch_fasta_translated    // channel : [ path(translated.fasta) ]
    ch_versions            // channel : [ path(versions.yml) ]


    main:

    // DESC: Query all of IDs from serialized JSON file and save them to a file
    // ARGS: qresults_serialized (channel)        - Channel containing path to serialized JSON file
    // RETS: ch_qresults_serialized_ids (channel) - Channel containing path to file containing all IDs
    //       ch_versions (channel)                - Channel containing path to versions.yml
    JQ_QUERY_ID(
        qresults_serialized = ch_qresults_serialized,
    )
    ch_versions = ch_versions.mix(JQ_QUERY_ID.out.versions)

    // DESC: Extract sequences from fasta file containing all transalted ORFs
    // ARGS: headers (channel)           - Channel containing path to file containing all IDs
    //       fasta (channel)             - Channel containing path to translated ORFs
    // RETS: ch_seqkit_grepseq (channel) - Channel containing path to file containing all sequences
    //       ch_versions (channel)       - Channel containing path to versions.yml
    SEQKIT_GREP(
        headers = JQ_QUERY_ID.out.query_ids,
        fasta = ch_fasta_translated,
    )
    ch_versions = ch_versions.mix(SEQKIT_GREP.out.versions)

    // DESC: Add sequences to the serialized JSON file
    // ARGS: qresults_serialized (channel)    - Channel containing path to serialized JSON file
    //       fasta (channel)                  - Channel containing path to file containing all sequences
    //       type (str)                       - Type of sequence to add
    // RETS: ch_qresults_serialized (channel) - Channel containing path to serialized JSON file with sequences
    //       ch_versions (channel)            - Channel containing path to versions.yml
    DOMTBLOP_ADDSEQ(
        qresults_serialized = ch_qresults_serialized,
        fasta = SEQKIT_GREP.out.seqkit_grepseq,
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
    // ARGS: ch_qresults_serialized (channel) - Channel containing path to serialized JSON file
    // RETS: ch_qresults_bed (channel)        - Channel containing path to bed file
    //       ch_versions (channel)            - Channel containing path to versions.yml
    DOMTBLOP_TOBED(
        qresults_serialized = ch_qresults_serialized,
    )
    ch_versions = ch_versions.mix(DOMTBLOP_TOBED.out.versions)


    // DESC: Get the sequences from the fasta file based on the bed file
    // ARGS: fasta (channel)                - Channel
    //       bed (channel)                  - Channel
    // RETS: ch_bedtools_getfasta (channel) - Channel
    //       ch_versions (channel)          - Channel
    BEDTOOLS_GETFASTA(
        fasta = ch_fasta,
        bed = DOMTBLOP_TOBED.out.qresults_bed,
    )
    ch_versions = ch_versions.mix(BEDTOOLS_GETFASTA.out.versions)


    // DESC: Add nucleotide sequences to the serialized JSON file
    // ARGS: qresults_serialized (channel) - Channel containing path to serialized JSON file
    //       fasta (channel)               - Channel containing path to file containing all sequences
    //       type (str)                    - Type of sequence to add
    // RETS: ch_qresults_serialized (channel) - Channel containing path to serialized JSON file with sequences
    //       ch_versions (channel)            - Channel containing path to versions.yml
    DOMTBLOP_ADDSEQ(
        qresults_serialized = ch_qresults_serialized,
        fasta = BEDTOOLS_GETFASTA.out.bedtools_getfasta,
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
    fasta                 // str     : path(path/to/nucleotide.fasta)
    ch_fasta_translated   // channel : [ path(path/to/nucleotide.translated.fasta) ]
    ch_hmmscan_domtblout  // channel : [ path(path/to/domtblout) ]
    domain_ievalue        // str     : E-value threshold for domain filtering
    fasta_type            // str     : Type of fasta file (dna, rna, protein)
    min_alignment_len     // int     : Minimum length of alignment
    group                 // int     : Hit maximum grouping distance
    ch_versions           // channel : [ path(versions.yml) ]

    main:

    // DESC: Parse and seralize the results of each query to JSON
    // ARGS: hmmscan_domtblout (channel)   -  Channel contianing path to hmmscan domtblout file
    // RETS: ch_qresults_serialized (channel) - Channel containing path to serialized JSON file
    //       ch_versions (channel)            - Channel containing path to versions.yml
    DOMTBLOP_PARSER(
        hmmscan_domtblout = ch_hmmscan_domtblout,
    )
    ch_versions = ch_versions.mix(DOMTBLOP_PARSER.out.versions)


    // DESC: Add sequence data to the serialized JSON file
    // ARGS: ch_qresults_serialized (channel) - Channel containing path to serialized JSON file
    //       ch_fasta_translated (channel)    - Channel containing path to translated ORFs
    // RETS: ch_qresults_serialized (channel) - Channel containing path to serialized JSON file with sequences
    //       ch_versions (channel)            - Channel containing path to versions.yml
    DOMTBLOP_ADD_AMINOACIDSEQ_WORKFLOW(
        ch_qresults_serialized = DOMTBLOP_PARSER.out.qresults_serialized,
        ch_fasta = ch_fasta_translated,
        ch_versions = ch_versions,
    )
    ch_versions = ch_versions.mix(DOMTBLOP_ADD_AMINOACIDSEQ_WORKFLOW.out.versions)


    if (fasta_type == "dna" || fasta_type == "rna") {
        // DESC: Add nucleotide sequences to the serialized JSON file
        // ARGS: ch_qresults_serialized (channel) - Channel containing path to serialized JSON file with sequences
        //       ch_fasta (channel)               - Channel containing path to nucleotide fasta file
        // RETS: ch_qresults_serialized (channel) - Channel containing path to serialized JSON file with sequences
        //       ch_versions (channel)            - Channel containing path to versions.yml
        DOMTBLOP_ADD_NUCLEOTIDESEQ_WORKFLOW(
            ch_qresults_serialized = DOMTBLOP_ADD_AMINOACIDSEQ_WORKFLOW.out.qresults_serialized,
            ch_fasta = fasta,
            ch_versions = ch_versions,
        )
        ch_versions = ch_versions.mix(DOMTBLOP_ADD_NUCLEOTIDESEQ_WORKFLOW.out.versions)
        ch_qresults_serialized = DOMTBLOP_ADD_NUCLEOTIDESEQ_WORKFLOW.out.qresults_serialized
    } else {
        ch_qresults_serialized = DOMTBLOP_ADD_AMINOACIDSEQ_WORKFLOW.out.qresults_serialized
    }


    // DESC: Filter hits using each domain individual E-value as a threshold
    // ARGS: ch_qresults_serialized (channel) - Channel containing path to serialized JSON file with sequences
    //       threshold (str)                 - E-value threshold (e.g. 1e-5, 0.001)
    // RETS: ch_qresults_serialized (channel) - Channel containing path to serialized JSON file with filtered hits
    //       ch_versions (channel)            - Channel containing path to versions.yml
    DOMTBLOP_FILTER_BY_DOMAIN_IEVALUE(
        ch_qresults_serialized = ch_qresults_serialized,
        threshold = domain_ievalue,
    )
    ch_versions = ch_versions.mix(DOMTBLOP_FILTER_BY_DOMAIN_IEVALUE.out.versions)


    // DESC: Filter hits by minimum alignment length
    // ARGS: ch_qresults_serialized (channel) - Channel containing path to serialized JSON file with sequences
    //       min_alignment_len (int)          - Minimum length of alignment
    // RETS: ch_qresults_serialized (channel) - Channel containing path to serialized JSON file with filtered hits
    //       ch_versions (channel)            - Channel containing path to versions.yml
    DOMTBLOP_FILTER_BY_MIN_ALIGNMENT_LENGTH(
        ch_qresults_serialized = DOMTBLOP_FILTER_BY_DOMAIN_IEVALUE.out.qresults_serialized,
        min_alignment_len,
    )
    ch_versions = ch_versions.mix(DOMTBLOP_FILTER_BY_MIN_ALIGNMENT_LENGTH.out.versions)


    // DESC: Group hits within specified distance from each other.
    // ARGS: ch_hits_serialized (channel) - Channel containing path to serialized JSON file with hits
    //       group (int)                  - Maximum distance between hits to group
    // RETS: ch_qresults_serialized (channel) - Channel containing path to serialized JSON file with grouped hits
    //       ch_versions (channel)            - Channel containing path to versions.yml
    DOMTBLOP_GROUP(
        ch_hits_serialized = DOMTBLOP_FILTER_BY_MIN_ALIGNMENT_LENGTH.out.qresults_serialized,
        group
    )
    ch_versions = ch_versions.mix(DOMTBLOP_GROUP.out.versions)


    // DESC: Convert domain aligments from serialized JSON to GFF format
    // ARGS: ch_qresults_serialized (channel) - Channel containing path to serialized JSON file
    // RETS: ch_gff (channel)                 - Channel containing path to GFF file
    //       ch_versions (channel)            - Channel containing path to versions.yml
    DOMTBLOP_TOGFF(
        DOMTBLOP_GROUP.out.qresults_serialized,
    )
    ch_versions = ch_versions.mix(DOMTBLOP_TOGFF.out.versions)


    // DESC: Convert domain aligments from serialized JSON to FASTA format
    //       Writes: i) FASTA file containing nucleotide sequences
    //               ii) FASTA file containing amino acid sequences
    //               iii) FASTA file containing the amino acid sequences that could be aligned against the HMM
    // ARGS: ch_qresults_serialized (channel) - Channel containing path to serialized JSON file
    // RETS: ch_fasta (channel)               - Channel containing path to FASTA files
    //       ch_versions (channel)            - Channel containing path to versions.yml
    DOMTBLOP_TOFASTA(
        DOMTBLOP_GROUP.out.qresults_serialized,
    )
    ch_versions = ch_versions.mix(DOMTBLOP_TOFASTA.out.versions)


    emit:
    qresults_serialized = DOMTBLOP_GROUP.out.qresults_serialized
    versions = ch_versions

}
