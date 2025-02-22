import groovy.yaml.YamlBuilder
include { EXSCAN } from './workflows/exscan'


help_message = """
exscan.nf - An extended sequence scanner
========================================
Explore key elements in a sequence using HMM profile scanning.

Usage:
    nextflow run exscan.nf --fasta <fasta> --hmmdb <hmm_db>
    nextflow run exscan.nf -params-file <yaml>

Required Arguments:
    --fasta <fasta>             : Input fasta file

    --hmmdb <hmmdb>             : Profile HMM database file
                                  Index files from hmmpress must also be located in the same directory

Alternative to Required Arguments:
    -params-file <yaml>         : YAML/JSON file with the parameters
                                  Mutually exclusive with --fasta and --hmmdb

Optional Arguments:
    --sequence_type <val>       : Type of sequences in the input `--fasta` file.
                                  Choices: {dna, protein} [defalut: ${params.sequence_type}]

    --dom_ieval_filter <val>    : Filter domain alignment by individual E-value [default: ${params.dom_ieval_filter}]
                                  Retains alignments with individual E-values **equal or less** than this value

    --keep-best-hit             : Retain only the best hit for each query sequence according to the independent E-value
                                  [default: ${params.keep_best_hit}]

    --min_alignment_len <val>   : Minimum length for domain alignments [default: ${params.min_alignment_len}]
                                  Retains alignments with length **equal or more** than this value

    --grouping_distance <val>   : Group query results within this many base pairs [default: ${params.grouping_distance}]

    --gff_intersect <gff>       : If provided, any feature in this GFF intersecting
                                  with a domain alignment will be retained.

    --keep_only_intersect       : Retain only query results that intersect with features
                                  of the provided GFF file.

    --outdir <outdir>           : Output directory [default: ${params.outdir}]

    --help                      : Print help message and exit

    --version                   : Print version and exit
"""

init_summary = """
E X S C A N   P I P E L I N E   v${params.manifest.version}
======================================
fasta               : ${params.fasta}
hmmdb               : ${params.hmmdb}
sequence_type       : ${params.sequence_type}
dom_ieval_filter    : ${params.dom_ieval_filter}
keep_best_hit       : ${params.keep_best_hit}
min_alignment_len   : ${params.min_alignment_len}
grouping_distance   : ${params.grouping_distance}
gff_intersect       : ${params.gff_intersect}
keep_only_intersect : ${params.keep_only_intersect}
outdir              : ${params.outdir}

--

Run as              : ${workflow.commandLine}
Started at          : ${workflow.start}
Config files        : ${workflow.configFiles}

--
"""
// container images : ${workflow.containerEngine}:${workflow.container}


// DESC: Validate input arguments and initialize pipeline, printing a small summary
// ARGS: None, uses variables defined at the beginning of the script
// OUTS: `$init_summary` as log message at `INFO` level
//       `$help_message` as stdout if `--help` flag is set
//       `$version` as stdout if `--version` flag is set
//       Proper error message and exit code if required arguments are missing
// RETS: None
def validateParams() {

    // `--help` and `--version` flags
    if (params.help) {
        println help_message
        System.exit(0)
    }
    if (params.version) {
        println "${params.manifest.name} v${params.manifest.version}"
        System.exit(0)
    }

    // Check required arguments
    if (params.fasta == null) {
        println help_message
        log.error "Missing required argument: --fasta"
        System.exit(1)
    }
    if (!file(params.fasta).exists()) {
        log.error "File not found: ${params.fasta}"
        System.exit(1)
    }

    if (params.hmmdb == null) {
        println help_message
        log.error "Missing required argument: --hmmdb"
        System.exit(1)
    }
    if (!file(params.hmmdb).exists()) {
        log.error "File not found: ${params.hmmdb}"
        System.exit(1)
    }

    if (params.sequence_type != 'dna' && params.sequence_type != 'protein') {
        log.error "Invalid fasta type: ${params.sequence_type}. Must be one of {dna, protein}"
        System.exit(1)
    }

}


// DESC: Handle hmmdb appropriately so it can be mounted correctly
// ARGS: None, uses variables defined at the beginning of the script
// OUTS: `$params.hmmdb_file` and `$params.hmmdb_dir` as new variables
// RETS: None
def handle_hmmdb() {
    // hmmscan needs access to the files generated by hmmpress,
    // but accepts the hmmdb file as input.
    // When using docker for example, if we only mount the hmmdb file,
    // hmmscan will not be able to find the index files.
    // To solve this, we mount the directory containing the hmmdb file
    // and pass the hmmdb file as input to hmmscan.
    params.hmmdb_file = params.hmmdb.split('/').last()
    params.hmmdb_dir = file(params.hmmdb).getParent()
}


// DESC: Display completion message based on workflow status
// ARGS: None, uses variables defined at the beginning of the script
// OUTS: Completion message at `INFO` or `ERROR` level
// RETS: None
def completionMsg() {

    if (workflow.success) {
        if (workflow.stats.ignoredCount == 0) {
            log.info "Pipeline completed successfully!"
        }
        else {
            log.info "Pipeline completed successully, but with errored processes"
        }
    }
    else {
        log.error "Pipeline completed with errors"
    }

}


// DESC: Dump parameters to a YAML file
// ARGS: None, uses variables defined at the beginning of the script
// OUTS: YAML file with the parameters
// RETS: None
def dumpParametersToYaml() {
    def paramMap = [
        fasta: params.fasta,
        hmmdb: params.hmmdb,
        sequence_type: params.sequence_type,
        dom_ieval_filter: params.dom_ieval_filter,
        keep_best_hit: params.keep_best_hit,
        min_alignment_len: params.min_alignment_len,
        grouping_distance: params.grouping_distance,
        gff_intersect: params.gff_intersect,
        keep_only_intersect: params.keep_only_intersect,
        outdir: params.outdir,
    ]

    def yaml = new YamlBuilder()
    yaml(paramMap)

    def outputDir = new File("${params.outdir}/pipeline_info")
    outputDir.mkdirs()
    def outputFile = new File(outputDir, "params_used.yaml")
    outputFile.text = yaml.toString()
}


// Main workflow
workflow {

    main:

    // Validate input parameters
    validateParams()
    // Handle hmm database paths
    handle_hmmdb()
    // Initialization Summary - Everything looks good so far
    log.info init_summary

    dumpParametersToYaml()


    ch_versions = Channel.empty()
    // WORKFLOW: After validation, main workflow is launched here
    EXSCAN(
        params.dom_ieval_filter,
        params.fasta,
        params.sequence_type,
        params.gff_intersect,
        params.grouping_distance,
        params.hmmdb_dir,
        params.hmmdb_file,
        params.keep_best_hit,
        params.keep_only_intersect,
        params.min_alignment_len,
        ch_versions,
    )
    ch_versions = ch_versions.mix(EXSCAN.out.versions)




    // Save versions of all tools used in the pipeline
    ch_versions.collectFile(
        storeDir: "${params.outdir}/pipeline_info/",
        name: 'versions.yml',
        sort: true,
        newLine: true
    )

    // Display any error encountered during the workflow
    workflow.onComplete {
        completionMsg()
    }
}
