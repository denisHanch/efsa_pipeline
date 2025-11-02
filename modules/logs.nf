
import java.time.format.DateTimeFormatter
import java.time.ZoneId

def logUnmapped(pct, threshold, out_folder_name) {
    pct.collect{
    def percentage = it.toDouble()
    def log_line = "The percentage of unmapped reads in ${out_folder_name} pipeline: ${percentage} %"

    if (percentage == 0.0)
        log.info log_line
    else if (percentage <= threshold)
        log.warn log_line
    else
        log.error log_line + "\nChange parameter map_to_mod_fa in nextflow.config to true."
    }
}

def describePipeline(read_type, fasta_type) {
    return "▶ Running pipeline processing ${read_type} reads - mapping to the ${fasta_type} fasta."
}


def logWorkflowCompletion(out_folder_name, run_flag) {
    workflow.onComplete {

    def formatter = DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss")
                                .withZone(ZoneId.systemDefault())
    def readableTime = formatter.format(workflow.complete)


        if (run_flag) {
            if (workflow.success) {
                log.info "✅ The ${out_folder_name} processing pipeline completed successfully. ${readableTime}\n"
            } else {
                log.error "❌ The ${out_folder_name} processing pipeline failed: ${workflow.errorReport}\n"
            }
        }
    }
}