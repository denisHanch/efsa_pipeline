
import java.time.format.DateTimeFormatter
import java.time.ZoneId

def logUnmapped(reads, total_reads, out_folder_name, reference) {
    reads.combine(total_reads).map { r, total ->
    def percentage = r.toInteger() * 100.0 / total.toInteger()
    def pctStr = String.format('%.2f', percentage)
    log.info "ℹ️ The number of unmapped reads in ${out_folder_name} pipeline ${reference}: ${r} (${pctStr} %). Total input: ${total})\n"
    }
}

def describePipeline(read_type, fasta_type) {
    return "▶ Running pipeline processing ${read_type} reads - mapping to the ${fasta_type} fasta."
}


def logWorkflowCompletion(out_folder_name) {
    workflow.onComplete {

        if (workflow.success) {

            log.info "✅ The ${out_folder_name} processing pipeline completed successfully."

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
                    log.info "ℹ️ Nextflow work/ directory was removed."
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