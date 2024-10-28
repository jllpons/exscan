#!/usr/bin/env nextflow

help_message = """
E X S C A N   P I P E L I N E
=============================
Exploration and annotation of exons and their features using profile HMMs.

Usage:
    nextflow run exscan.nf --fasta <fasta> --hmm_db <hmm_db>

Required arguments:
    --fasta <fasta>     : Input fasta file
    --hmmdb <hmmdb>     : Profile HMM database directory, where the hmmpress
                          index files are located.
                          **The file that was passed to hmmpress must be
                          also located in the same directory as the index files
                          and must be named 'hmmdb'**

Optional arguments:
    --outdir <outdir>   : Output directory [default: results]
    --help              : Print help message and exit
    --version           : Print version and exit
"""

version = "0.0.0"

params.fasta = null
params.hmmdb = null
params.outdir = 'results'

params.help = false
params.version = false

init_summary = """
E X S C A N   P I P E L I N E
=============================
fasta            : ${params.fasta}
profile database : ${params.hmmdb}
outdir           : ${params.outdir}
"""

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { WRITE_VERSIONS } from './modules/local/utils'
include { EXSCAN } from './workflows/exscan'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// DESC: Validate input arguments and initialize pipeline
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
        println "exscan.nf ${version}"
        System.exit(0)
    }

    // Check required arguments
    if (params.fasta == null) {
        println help_message
        log.error "Missing required argument: --fasta"
        System.exit(1)
    }
    if (!file(params.fasta).exists()) {
        log.error "File not found: ${fasta}"
        System.exit(1)
    }

    if (params.hmmdb == null) {
        println help_message
        log.error "Missing required argument: --hmmdb"
        System.exit(1)
    }
    if (!file(params.hmmdb).exists()) {
        log.error "File not found: ${hmmdb}"
        System.exit(1)
    }

    // Everything looks good so far
    log.info init_summary
}

def completionTasks() {

    if (workflow.success) {
        if (workflow.stats.ignoredCount == 0) {
            log.info """
Pipeline completed successfully
            """
        }
        else {
            log.info """
Pipeline completed successully, but with erroed processes
            """
        }
    }
    else {
        log.error """
Pipeline completed with errors
        """
    }

}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/



workflow {

    main:

    validateParams()

    // Channel for input fasta file
    ch_fasta = Channel.fromPath(params.fasta)
    // Channel for profile HMM database directory
    ch_hmmdb = Channel.fromPath(params.hmmdb)
    // Channel where each tool will dump its version information
    ch_versions = Channel.empty()


    // WORKFLOW: After validation, main workflow is launched here
    EXSCAN(
        ch_fasta,
        ch_hmmdb,
        ch_versions
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
        completionTasks()
    }
}
