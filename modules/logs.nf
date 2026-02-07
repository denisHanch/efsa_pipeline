
import java.time.format.DateTimeFormatter
import java.time.ZoneId

def logUnmapped(reads, total_reads, out_folder_name, reference) {

    reads.combine(total_reads).subscribe { r, total ->

        long unmapped = r as long
        long totalInput = total as long

        def percentage = (unmapped * 100.0) / totalInput
        def pctStr = String.format("%.2f", percentage)

        // Build a plain String
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

        if (workflow.success) {
            def formatter = DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss").withZone(ZoneId.systemDefault())
            def readableTime = formatter.format(workflow.complete)

            log.info "✅ The ${out_folder_name} processing pipeline completed successfully.\n"

            def workDir = new File("${workflow.workDir}")
            def launchDir = new File("${workflow.launchDir}")
            def logDir = new File("${params.log_dir}")
            logDir.mkdirs()

            workDir.eachFileRecurse(groovy.io.FileType.ANY) { f ->
            if( f.name ==~ /^(\.command).*/ ) {
                
                def relPath = workDir.toPath().relativize(f.toPath()).toString()
                def dest = new File(logDir, relPath)
                dest.parentFile.mkdirs()
                f.withInputStream { ins -> dest.withOutputStream { out -> out << ins } }
                }
            }

            if (params.clean_work && out_folder_name == "execution of main.nf") {
                if( workDir.exists() ) {
                    workDir.deleteDir()
                    log.info "ℹ️  Nextflow work/ directory was removed.\n"
                }
            }

        } else {
            log.error "❌ The ${out_folder_name} processing pipeline failed: ${workflow.errorReport}"
        }
    }
}


def loadFastqFiles(pathPattern) { 
    return Channel.fromPath(pathPattern)
                  .map { file ->
                      def name = file.baseName.replaceFirst(/\.fastq/, '')
                      [name, file]
                  }
}

def loadShortFastqFiles(short_read_files) {
    return Channel.from(short_read_files).map { file ->
            def matcher = file.name =~ /^(.+?)(?:[_\.](S[0-9]+_L[0-9]+_)?(R[12]|[12]))?\.f(ast)?q\.gz$/
            if (matcher.matches()) {
                [matcher[0][1], file]
            }
        }
        .groupTuple(sort: true)
}

def listFastqFiles(String dirPath, String pattern = ".*\\.(fastq|fq)(\\.gz)?\$") {
    def dir = file(dirPath)
    if (!dir.exists() || !dir.isDirectory()) {
        log.warn "Directory not found or not a directory: ${dirPath}"
        return []
    }
    def files = dir.listFiles()?.findAll { it.name =~ pattern } ?: []
    if (!files) log.info "No FASTQ files found in ${dirPath} matching pattern ${pattern}"
    return files
}