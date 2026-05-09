#!/usr/bin/env nextflow

include { validate } from "./modules/validate.nf"
include { analysis } from "./workflows/analysis.nf"
include { generateProcessManifest ; logToNextflowFile ; logCompletionSummary ; getLogDir ; copyCommandLogs ; resolveWorkDir } from "./modules/logs.nf"

// Help message
def helpMessage() {
    log.info(
        """
    Usage:

    nextflow run main.nf
    
    Options:

    -resume          Run pipeline from the point where it was interrupted or failed (Nextflow built-in)
    --config_json    Path to the JSON file containing parameters for the pipeline   (default: ${params.config_json})
    --out_dir        Output directory                                               (default: ${params.out_dir})
    --max_cpu        Maximum CPUs per process                                       (default: ${params.max_cpu})
    --clean_work     Remove workdir after success                                   (default: ${params.clean_work})
    --validation-level <level> Validation strictness: STRICT, TRUST, or MINIMAL     (default: ${params.validation_level}, overridden by config.json)
    --logging-level <level>    Log verbosity: DEBUG, INFO, WARNING, or ERROR        (default: ${params.logging_level}, overridden by config.json)
    --organism_type <type>     Organism type: PROKARYOTE or EUKARYOTE               (default: ${params.organism_type}, overridden by config.json)
    --force-defragment-ref     Force reference defragmentation [UNSUPPORTED]        (overridden by config.json)
    -with-report     Generate HTML execution report                                 (Nextflow built-in)
    -with-timeline   Produce timeline visualization                                 (Nextflow built-in)
    -with-dag        Produce DAG of workflow                                        (Nextflow built-in)
    --help           Show this help message
    """.stripIndent()
    )
}


workflow {
    // Show help
    if (params.help) {
        helpMessage()
        exit(0)
    }

    file("${params.out_dir}/tables/csv_per_sv_summary").mkdirs()

    config_ch = channel.fromPath(params.config_json, checkIfExists: true)
    validate(config_ch)
    analysis(validate.out.params_json)

    workflow.onComplete { wf ->
        def workDir = resolveWorkDir(wf)
        def logDir = getLogDir()

        copyCommandLogs(workDir, logDir)
        generateProcessManifest(logDir, wf)
        logCompletionSummary(wf, workDir)

        if (workDir.deleteDir()) {
            logToNextflowFile("🧹 Removed work directory: ${workDir.absolutePath}\n")
        }
        else {
            logToNextflowFile("⚠️ Failed to remove work directory: ${workDir.absolutePath}\n")
        }
    }

    workflow.onError {
        def errorDetails = workflow?.errorMessage ?: workflow?.errorReport ?: 'No error details available (pipeline may have been interrupted)'
        def logDirPath = (params?.log_dir ?: 'data/outputs/logs')
        logToNextflowFile("Pipeline execution stopped with the following message: ${errorDetails}")
        logToNextflowFile("Check the process execution manifest in ${logDirPath}/process_manifest.txt for details on which processes failed.")
    }
}
