
import java.time.format.DateTimeFormatter
import java.time.ZoneId

def logUnmapped(reads, total_reads, out_folder_name, reference) {

    reads.combine(total_reads).subscribe { r, total ->

        long unmapped = r as long
        long totalInput = total as long

        def percentage = (unmapped * 100.0) / totalInput
        def pctStr = String.format("%.2f", percentage)

        String msg = "📊 ${out_folder_name} mapping${reference}:\n" +
                     "    Unmapped reads: ${String.format("%,d", unmapped)} (${pctStr}%)\n" +
                     "    Total input reads: ${String.format("%,d", totalInput)}\n"

        log.info(msg)
    }
}


def describePipeline(read_type, fasta_type) {
    log.info "ℹ️  Running pipeline: processing ${read_type} reads → mapping to the ${fasta_type} fasta.\n"
}



def logWorkflowCompletion(out_folder_name) {
    workflow.onComplete {
        def formatter = DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss").withZone(ZoneId.systemDefault())
        def readableTime = formatter.format(workflow.complete)

        def workDir = new File("${workflow.workDir}")
        def launchDir = new File("${workflow.launchDir}")
        def logDir = new File("${params.log_dir}")
        logDir.mkdirs()

        // Always copy process logs to output directory (both on success and failure)
        if (workDir.exists()) {
            workDir.eachFileRecurse(groovy.io.FileType.ANY) { f ->
                if( f.name ==~ /^(\.command).*/ ) {
                    def relPath = workDir.toPath().relativize(f.toPath()).toString()
                    def dest = new File(logDir, relPath)
                    dest.parentFile.mkdirs()
                    f.withInputStream { ins -> dest.withOutputStream { out -> out << ins } }
                }
            }
        }

        // Generate process execution manifest with success/failure status and exit codes
        generateProcessManifest(logDir)

        if (workflow.success) {
            log.info "✅ The ${out_folder_name} processing pipeline completed successfully.\n"

            if (params.clean_work && out_folder_name == "execution of main.nf") {
                if ( workDir.exists() ) {
                    workDir.deleteDir()
                    log.info "ℹ️  Nextflow work/ directory was removed.\n"
                }
            }
        } else {
            log.error "❌ The ${out_folder_name} processing pipeline failed: ${workflow.errorReport}"
        }
    }
}

/**
 * Parse the Nextflow trace file and produce a human-readable manifest
 * listing every process execution with its status and exit code.
 */
def generateProcessManifest(File logDir) {
    def traceFile = new File(logDir, "trace.tsv")
    def manifestFile = new File(logDir, "process_manifest.txt")

    manifestFile.text  = "# Pipeline Execution Manifest\n"
    manifestFile.append("# Generated: ${new Date()}\n")
    manifestFile.append("# Pipeline status: ${workflow.success ? 'SUCCESS' : 'FAILED'}\n")
    manifestFile.append("# Duration: ${workflow.duration}\n")
    if (workflow.errorMessage) {
        manifestFile.append("# Error: ${workflow.errorMessage}\n")
    }
    manifestFile.append("#\n")

    if (traceFile.exists()) {
        def lines = traceFile.readLines()
        // Copy full trace data into the manifest
        lines.each { line -> manifestFile.append(line + '\n') }

        if (lines.size() > 1) {
            def dataLines = lines.drop(1)
            def completed = dataLines.count { it.split('\t')[2]?.trim() == 'COMPLETED' }
            def failed    = dataLines.count { it.split('\t')[2]?.trim() == 'FAILED' }
            def cached    = dataLines.count { it.split('\t')[2]?.trim() == 'CACHED' }
            def aborted   = dataLines.size() - completed - failed - cached

            manifestFile.append("\n# ── Summary ──\n")
            manifestFile.append("# Completed: ${completed}\n")
            manifestFile.append("# Failed:    ${failed}\n")
            manifestFile.append("# Cached:    ${cached}\n")
            manifestFile.append("# Aborted:   ${aborted}\n")
            manifestFile.append("# Total:     ${dataLines.size()}\n")

            def failedProcesses = dataLines.findAll { it.split('\t')[2]?.trim() == 'FAILED' }
            if (failedProcesses) {
                manifestFile.append("\n# ── Failed Processes ──\n")
                failedProcesses.each { line ->
                    def fields = line.split('\t')
                    manifestFile.append("# Process: ${fields[1]}, Exit code: ${fields[3]}\n")
                }
            }
        }
    } else {
        manifestFile.append("# Note: Trace file not found. Enable trace in nextflow.config for process-level details.\n")
    }

    log.info "📋 Process execution manifest: ${manifestFile.path}\n"
}