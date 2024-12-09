
include { GAWK_HEADER_METADATA    } from '../modules/local/gawk'
include { HMMSCAN                 } from '../modules/local/hmmscan'
include { JQ_QUERY_ID             } from '../modules/local/jq'
include { MERGE_DOMTBLOUT_RESULTS } from '../modules/local/hmmscan'
include { SEQKIT_GREP             } from '../modules/local/seqkit'
include { SEQKIT_SPLIT            } from '../modules/local/seqkit'
include { SEQKIT_TRANSLATE        } from '../modules/local/seqkit'

include { DOMTBLOP_DEFAULT } from '../subworkflows/domtblop'

workflow EXSCAN {

    take:
    domain_ieval        // (float) - evalue threshold for domain hits
    fasta               // (str) - path to nucleotide fasta file
    fasta_type          // (str) - type of fasta file (dna, rna, or protein)
    gff_intersect       // (str) - path to gff file to intersect with hmmscan results
    group_distance      // (int) - distance in bp under which two hmmscan hits are grouped
    hmmdb_dir           // (str) - path to HMM database directory
    hmmdb_file          // (str) - Name of HMM database file
    keep_only_intersect // (bool) - keep only hmmscan hits that intersect with features in gff_intersect
    min_alignment_len   // (int) - minimum length of alignment
    ch_versions         // (channel) - channel containing path to versions.yml

    main:

    ch_fasta = Channel.fromPath(fasta)


    if (fasta_type == 'dna') {
        // DESC: Translate all 6 ORFs of a nucleotide fasta file.
        //       A fasta file containing the translated ORFs is returned.
        //       One fasta entry = whatever translated sequence bewteen two stop codons.
        // ARGS: fasta (str)                   - path to nucleotide fasta file
        // RETS: ch_fasta_translated (channel) - channel containing path to translated ORFs
        //       ch_versions (channel)         - channel containing path to versions.yml

        SEQKIT_TRANSLATE(
            fasta      = ch_fasta,
            fasta_type = fasta_type
        )
        // After each process, software versions are collected
        ch_versions = ch_versions.mix(SEQKIT_TRANSLATE.out.versions)
        ch_fasta_translated = SEQKIT_TRANSLATE.out.fasta_translated
    } else if (fasta_type == 'rna') {
        SEQKIT_TRANSLATE(
            fasta      = ch_fasta,
            fasta_type = fasta_type
        )
        ch_versions = ch_versions.mix(SEQKIT_TRANSLATE.out.versions)
        // DESC: If the fasta file is already translated, then we can skip the translation step.
        //       Still, we need to mock the metadata that seqkit-translate would have appended
        //       to the fasta headers.
        // ARGS: fasta (str)                   - path to nucleotide fasta file
        // RETS: ch_fasta_translated (channel) - channel containing path to translated ORFs
        //       ch_versions (channel)         - channel containing path to versions.yml
        GAWK_HEADER_METADATA(
            fasta = SEQKIT_TRANSLATE.out.fasta_translated
        )
        ch_versions = ch_versions.mix(GAWK_HEADER_METADATA.out.versions)
        ch_fasta_translated = GAWK_HEADER_METADATA.out.fasta_translated
    } else {
        GAWK_HEADER_METADATA(
            fasta = ch_fasta
        )
        ch_versions = ch_versions.mix(GAWK_HEADER_METADATA.out.versions)
        ch_fasta_translated = GAWK_HEADER_METADATA.out.fasta_translated
    }


    // DESC: Split large fasta files into smaller files so that they can be processed in parallel by hmmscan.
    //       Note: hmmscan (as far as I know) is very I/O intensive, so parallel > multithreaded.
    // ARGS: fasta (str)                   - path to nucleotide fasta file
    // RETS: ch_fasta_translated (channel) - channel containing path to translated ORFs
    //       ch_versions (channel)         - channel containing path to versions.yml
    SEQKIT_SPLIT(
        fasta = ch_fasta_translated
    )
    ch_versions = ch_versions.mix(SEQKIT_SPLIT.out.versions)


    // DESC: Query each ORF against a domain profile HMM database to identify
    //       potential homologous domains.
    // ARGS: hmmdb_file (str)               - Name of HMM database file
    //       hmmdb_dir (str)                - path to HMM database directory
    //       ch_fasta_translated (channel)  - channel containing path to a structured ORF file
    // RETS: ch_hmmscan_domtblout (channel) - path to hmmscan domtblout file
    HMMSCAN(
        hmmdb_file  = hmmdb_file,
        hmmdb_dir   = hmmdb_dir,
        fasta       = SEQKIT_SPLIT.out.fasta_split.flatten()
    )
    ch_versions = ch_versions.mix(HMMSCAN.out.versions)


    Channel.empty()
        .mix( HMMSCAN.out.hmmscan_domtblout )
        .collect()
        .set { ch_hmmscan_domtblout_results }
    // DESC: Large fasta_translated files are split into multiple files,
    //       so hmmscan can be run in parallel.
    //       This process will merge the results back into a single file.
    // ARGS: HMMSCAN.out.hmmscan_domtblout (channel) - channel containing path to hmmscan domtblout file
    // RETS: ch_hmmscan_merged_domtblout (channel)   - channel containing path to merged hmmscan domtblout file
    MERGE_DOMTBLOUT_RESULTS(
        ch_hmmscan_domtblout = ch_hmmscan_domtblout_results
    )


    // DESC: Parse hmmscan domtblout and perform different opeations. (See `./subworkflows/domtblop.nf`)
    // ARGS: fasta (str)                       - path to nucleotide fasta file
    //       ch_fasta_translated (channel)     - channel containing path to translated ORFs
    //       ch_hmmscan_domtblout (channel)    - cannel containing path to hmmscan domtblout file
    //       domain_ievalue (float)            - evalue threshold for domain hits
    //       fasta_type (str)                  - type of fasta file (dna, rna, or protein)
    //       gff_intersect (str)               - path to gff file to intersect with hmmscan results
    //       keep_only_intersect (bool)        - keep only hmmscan hits that intersect with gff file(s)
    //       min_alignment_len (int)           - minimum length of alignment
    //       group (int)                       - distance in bp under which two hmmscan hits are grouped
    //       ch_versions (channel)             - channel containing path to versions.yml
    // RETS: ch_qresults_serialized (channel)  - channel containing path to domtblop results
    DOMTBLOP_DEFAULT(
        fasta                = ch_fasta,
        ch_fasta_translated  = ch_fasta_translated,
        ch_hmmscan_domtblout = MERGE_DOMTBLOUT_RESULTS.out.hmmscan_merged_domtblout,
        domain_ieval         = domain_ieval,
        fasta_type           = fasta_type,
        gff_intersect        = gff_intersect,
        keep_only_intersect  = keep_only_intersect,
        min_alignment_len    = min_alignment_len,
        group_distance       = group_distance,
        ch_versions          = ch_versions
    )
    ch_versions = ch_versions.mix(DOMTBLOP_DEFAULT.out.versions)


    emit:
    results  = DOMTBLOP_DEFAULT.out.qresults_serialized // channel containing path to domtblop results
    versions = ch_versions                              // channel containing path to versions.yml
}

