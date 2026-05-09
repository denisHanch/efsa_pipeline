
def logToNextflowFile(String message) {
    def logDirPath = (params?.log_dir ?: 'data/outputs/logs').toString()
    def logDir = new File(logDirPath)
    logDir.mkdirs()
    def nfLog = new File(logDir, "nextflow.log")
    nfLog << "${message}${message.endsWith('\n') ? '' : '\n'}"
}

def logUnmapped(reads, total_reads, out_folder_name, reference) {
    reads.combine(total_reads).subscribe { r, total ->
        long unmapped = r as long
        long totalInput = total as long

        def percentage = (unmapped * 100.0) / totalInput
        def pctStr = String.format("%.2f", percentage)

        String msg = "📊 ${out_folder_name} mapping${reference}:\n" +
                     "    Unmapped reads: ${String.format("%,d", unmapped)} (${pctStr}%)\n" +
                     "    Total input reads: ${String.format("%,d", totalInput)}\n"

        logToNextflowFile(msg)
    }
}

/**
 * Parse the Nextflow trace file and produce a human-readable manifest
 * listing every process execution with its status and exit code.
 */
def generateProcessManifest(File logDir, wfMeta) {
    def traceFile = new File(logDir, "trace.tsv")
    def manifestFile = new File(logDir, "process_manifest.txt")

    manifestFile.text  = "# Pipeline Execution Manifest\n"
    manifestFile.append("# Generated: ${new Date()}\n")
    manifestFile.append("# Pipeline status: ${wfMeta?.success ? 'SUCCESS' : 'FAILED'}\n")
    manifestFile.append("# Duration: ${wfMeta?.duration ?: 'N/A'}\n")
    if (wfMeta?.errorMessage) {
        manifestFile.append("# Error: ${wfMeta.errorMessage}\n")
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

    logToNextflowFile("📋 Process execution manifest: ${manifestFile.path}\n")
}

def getLogDir() {
    def logDirPath = (params?.log_dir ?: 'data/outputs/logs').toString()
    def logDir = new File(logDirPath)
    logDir.mkdirs()
    return logDir
}

def resolveWorkDir(wfMeta) {
    def wfWorkDir = wfMeta?.workDir?.toString()
    if (wfWorkDir) {
        return new File(wfWorkDir)
    }

    def launchDirPath = wfMeta?.launchDir?.toString() ?: System.getProperty('user.dir')
    return new File(launchDirPath, "work")
}

def copyCommandLogs(File workDir, File logDir) {
    if (!workDir?.exists()) return

    workDir.eachFileRecurse(groovy.io.FileType.ANY) { f ->
        if (f.name ==~ /^(\.command).*/) {
            def relPath = workDir.toPath().relativize(f.toPath()).toString()
            def dest = new File(logDir, relPath)
            dest.parentFile.mkdirs()
            f.withInputStream { ins -> dest.withOutputStream { out -> out << ins } }
        }
    }
}

def logCompletionSummary(wf, File workDir) {
    if (wf?.success) {
        logToNextflowFile("✅ The execution of main.nf processing pipeline completed successfully.\n")
    }
}
